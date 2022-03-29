//
// Created by Tom Koller on 09.07.20.
//

#ifndef ZAVI_BOULDERN_BATCHESTIMATOR_H
#define ZAVI_BOULDERN_BATCHESTIMATOR_H
#include <ceres/ceres.h>
#include "ADEKF/ManifoldCreator.h"
#include "ADEKF/types/SO3.h"
#include "Plugins/FileIMUPlugin.hpp"
#include "Plugins/ARTPlugin.hpp"
#include "Plugins/SynchronizedSensorsPlugin.hpp"
//#include "../../BatchViterbi/include/ladder_functions.h"
#include <ADEKF/viz/adekf_viz.h>
#include <filesystem>
#include "oriestimu_src/CsgOriEstIMU.hpp"

/**
 * @brief Contains functions for the transition estimator
 * 
 */
namespace velocity_batch
{
    ADEKF_MANIFOLD(Velo3D, ((adekf::SO3, orientation)), (3, velocity))
    
    /**
     * @brief The template state type
     * 
     */
    template <typename T>
    using State = Velo3D<T>;

    /**
     * @brief The double state type 
     * 
     */
    using STATE_TYPE = State<double>;

    
    /**
     * @brief The dof of the state
     * 
     */
    constexpr size_t DOF = adekf::DOFOf<STATE_TYPE>;
    /**
     * @brief The number of parameters in the state
     * 
     */
    constexpr size_t STATE_SIZE = STATE_TYPE::GLOBAL_SIZE;
   
    /**
     * @brief Shorthand for double vectors
     * 
     * @tparam size size of the vectors
     */
    template <int size>
    using VecD = Eigen::Matrix<double, size, 1>;

/**
 * @brief Maps a vector from an array given the start position FROM
 * 
 */
#define MAP_VECTOR_FROM(ptr, size, TYPE, CONST, FROM) Eigen::Map<CONST Eigen::Matrix<TYPE, size, 1>>(&(ptr)[FROM])
/**
 * @brief Maps a vector from a pointer
 * 
 */
#define MAP_VECTOR(ptr, size, TYPE, CONST) MAP_VECTOR_FROM(ptr, size, TYPE, CONST, 0)
/**
 * @brief Retrieves only the velocity from a  state pointer
 * 
 */
#define getVel(name, TYPE, CONST) MAP_VECTOR_FROM(name, 3, TYPE, CONST, 4)
/**
 * @brief Retrieves the accelerometer bias
 * 
 */
#define getAccBias(name, TYPE, CONST) MAP_VECTOR_FROM(name, 3, TYPE, CONST, 0)
/**
 * @brief Retrieves the gyrometer bias
 * 
 */
#define getArBias(name, TYPE, CONST) MAP_VECTOR_FROM(name, 3, TYPE, CONST, 3)
/**
 * @brief retrieves the orientation from a state pointer
 * 
 */
#define getOrient(name) \
    adekf::SO3 { name }

/**
 * @brief The dynamic model of the INS system
 * 
 * @tparam T Scalar type for auto differentiation
 * @param state The state x_k
 * @param bias The gyrometer and accelerometer bias
 * @param body_angular_rate  the measured angular velocity in body coordinates
 * @param body_acceleration  the measured accelerationin body coordinates
 * @param time_diff the time difference since the last state
 * @return State<T> the propagated state x_k+1
 */
    template <typename T>
    State<T> dynamicModel(const T *const state, const T *const bias, const VecD<3> &body_angular_rate, const VecD<3> &body_acceleration, double time_diff)
    {
        State<T> result;
        adekf::SO3<T> orient = getOrient(state) * adekf::SO3{(body_angular_rate - getArBias(bias, T, const)) * time_diff};
        result.orientation = orient;
        Eigen::Matrix<T, 3, 1> world_acceleration = (orient * (body_acceleration - getAccBias(bias, T, const))) - zavi::plugin::IMU_Plugin::gravity;
        result.velocity = getVel(state, T, const) + world_acceleration * time_diff;
        return result;
    }

   /**
    * @brief Cost function for the dynamic model
    * 
    */
    struct BoulderDynamic
    {
        const VecD<3> body_angular_rate, body_acceleration;
        const double t_step;
        //Diagonal of the stiffness matrix
        Eigen::Matrix<double, DOF, 1> quotient_matrix;

        BoulderDynamic(const VecD<3> &body_angular_rate, const VecD<3> &body_acceleration, double t_step, Eigen::Matrix<double,DOF,1> stds) : body_angular_rate(body_angular_rate),
                                                                                                            body_acceleration(body_acceleration), t_step(t_step), quotient_matrix(stds)
        {

            assert(t_step > 0.);
        }
        /**
         * @brief Calculates the dynamic cost for ceres
         * 
         * @tparam T Scalar type for auto differentiation
         * @param past the state pointer of x_k
         * @param next the state pointer of x_k+1
         * @param bias the bias pointer
         * @param residual the dynamic cost (the output)
         * @return true always
         */
        template <typename T>
        bool operator()(const T *const past, const T *const next, const T *const bias, T *residual) const
        {
            (MAP_VECTOR(residual,DOF, T,)) = quotient_matrix.cwiseProduct(State<T>(next) - dynamicModel<T>(past, bias, body_angular_rate, body_acceleration, t_step));
            return true;
        }
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };

    /**
     * @brief Cost function for zero velocity
     * 
     */
    struct ZeroVelocityCost
    {
        //stiffness of the ZUPT
        static inline double stiffness = 1.;

        ZeroVelocityCost() {}

        /**
         * @brief Operator to calculate the cost
         * 
         * @tparam T Scalar type for auto differentiation
         * @param state the state pointer of x_k
         * @param residual the residual cost
         * @return true always
         */
        template <typename T>
        bool operator()(const T *const state, T *residual) const
        {

            Eigen::Matrix<T, 3, 1> diff = getVel(state, T, const);; // target is 0 vector
            (MAP_VECTOR(residual, 3, T, )) = diff * stiffness;         // TODO sigma scaling

            return true;
        }
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };

    /**
     * @brief Cost for the magnetometer measurement
     * 
     */
    struct MagnetometerCost
    {
        //magnetometer measurement at time k
        const VecD<3> mag_meas;
        //World magnetic field direction
        static inline VecD<3> mag_direction;
        //Stiffness of the cost function
        static inline double stiffness = 1.;

        MagnetometerCost(const VecD<3> &mag_meas) : mag_meas(mag_meas)
        {
        }
        /**
         * @brief Calculates the cost of the magnetometer measurement.
         * Corrects the measurement by the magnetometer bias and normalizes it. 
         * The corrected measurement is compared with the world magnet direction.
         * 
         * @tparam T Scalar type for auto differentiation
         * @param state The state pointer x_k
         * @param magnetometer_bias The pointer of the magnetometer bias
         * @param residual The cost (output)
         * @return true always
         */
        template <typename T>
        bool operator()(const T *const state, const T* const magnetometer_bias, T *residual) const
        {
            // Q*(meas-bias)-world_magnet
            Eigen::Matrix<T, 3, 1> corrected_meas = (mag_meas - MAP_VECTOR(magnetometer_bias, 3, T, const)).normalized();
            //Eigen::Matrix<double, 3, 1> corrected_meas=mag_meas.normalized();
            Eigen::Matrix<T, 3, 1> diff = getOrient(state) * corrected_meas - mag_direction;
            (MAP_VECTOR(residual, 3, T, )) = diff * stiffness; // TODO sigma scaling

            return true;
        }
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };

    /**
     * @brief Prior on the magnetometer bias to prevent unobservability
     * 
     */
    struct BiasCost{
        //Stiffness of the bias
         static inline double stiffness = 1.;
        template <typename T>
        bool operator()(const T* const magnetometer_bias, T *residual) const
        {
            (MAP_VECTOR(residual, 3, T, )) = MAP_VECTOR(magnetometer_bias,3,T,const) * stiffness; // TODO sigma scaling

            return true;
        }
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };

    /**
     * @brief Computes the transition estimate for a single hand. 
     * 
     * @param imu The IMU Plugin for the imu data
     * @param start_time The start time of the trial 
     * @param end_time The end time of the trial
     * @param is_grip  A boolean vector which determines when a grip happens
     * @param magnet_in_art The world magnet vector
     * @return a vector with all velocities, covariances and orientations of the transition sampling 
     */
    adekf::aligned_vector<std::tuple<VecD<3>,Eigen::Matrix<double,6,6>,adekf::SO3d>> transitionEstimation(std::shared_ptr<zavi::plugin::FileIMUPlugin> imu, double start_time, double end_time, std::vector<bool> &is_grip, const VecD<3> &magnet_in_art)
    {
        //Read stiffness from config file
        ZeroVelocityCost::stiffness = 1./(double)cfg->lookup("velocity_std"); 
        MagnetometerCost::stiffness = 1. / (double)cfg->lookup("magnet_std");
        BiasCost::stiffness = 1. / (double)cfg->lookup("magnet_bias_std");
        // Build the problem.
        ceres::Problem problem;

        //Configure problem
        ceres::Solver::Options options;
        options.minimizer_type = ceres::TRUST_REGION;
        options.trust_region_strategy_type = ceres::DOGLEG;
        options.preconditioner_type = ceres::CLUSTER_JACOBI;
        options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
        options.dogleg_type = ceres::SUBSPACE_DOGLEG;
        options.minimizer_progress_to_stdout = false;
        options.num_threads = 8;
        options.gradient_tolerance = 1e-10;
        options.max_num_iterations = cfg->lookup("velocity_estimation_iterations");
        // options.update_state_every_iteration=true;
        options.use_nonmonotonic_steps = true;
        options.check_gradients = false;
        //Set world magnet direction
        MagnetometerCost::mag_direction = VecD<3>{ magnet_in_art};
        LOG(INFO) << "Magnet direction: " << MagnetometerCost::mag_direction.transpose() LOG_END;

        LOG(INFO) << "Load all IMU data" LOG_END;
        imu->setTimeRange(start_time, end_time);
        LOG(INFO)<< "Time Range: " << start_time << " - " << end_time << std::endl;
        //allocate memory for all states
        size_t rowCount = (end_time - start_time) * (double)cfg->lookup("imu_freq");
        auto x = new double[STATE_SIZE * ((int)rowCount + 7)];
       
        //LOG(INFO) << "Setting standard states. Rowcount is: " << rowCount LOG_END;
        if (x==NULL)
            LOG(FATAL) << "Could not allocate memory. " LOG_END;
        //Set start states
        memset(x, 0, STATE_SIZE*sizeof(double));

        x[3] = 1;
        
        // time index
        size_t t = 0;
        
        //vector to store all times
        std::vector<double> times(rowCount + 7);
        times[0] = imu->nextTime;
        //bias pointer for accelerometer and angular rate
        double bias[6];
        //set to 0
        memset(bias, 0, sizeof(double) * 6);
        //bias pointer of magnetometer
        double magnetometer_bias[3];
        //set to 0
        memset(magnetometer_bias, 0, 3 * sizeof(double));
        //start orientation for CCsgOriEst Orientation estimator
        double init_q[4]{1,0,0,0};
        //orientation estimator for initial guess
        CCsgOriEstIMU orient_guess{cfg->lookup("imu_freq"),init_q,magnetometer_bias,true};
        //Set magnet direction for orientation estimator
        memcpy(&orient_guess.magrefImuframe[1],magnet_in_art.data(),sizeof(double)*3);
        //Acceleration storage for visualization
        adekf::aligned_vector<Eigen::Vector3d> accelerations;
        /**
         * @brief Callback to read all data from the IMU and to add cost functions to the ceres problem
         * 
         */
        auto add_costs = [&accelerations,&orient_guess,&t, &x, &problem, &rowCount, &times, &bias, &magnetometer_bias, &is_grip](zavi::plugin::FileIMUPlugin *imu, void *, double time_diff)
        {
            times[t + 1] = imu->nextTime;

            
            accelerations.push_back(imu->getMagnetometer());

            //Current stds for stiffness matrix
            Eigen::Matrix<double,DOF,1> stds;
            stds << (imu->getAngularVelocityCov().diagonal().cwiseSqrt()  +Eigen::Vector3d::Ones()*cfg->lookup("ar_nonlinearity")*imu->getAngularVelocity().cwiseAbs().maxCoeff() )* time_diff,
            (imu->getAccelerationCov().diagonal().cwiseSqrt() +Eigen::Vector3d::Ones()*cfg->lookup("acc_nonlinearity")*imu->getAcceleration().cwiseAbs().maxCoeff())* time_diff;
            //Error Handling for out of bounds IMU measurements. Increases covariance extremely
            if (imu->reading_error)
            {
                stds.segment<3>(3)+=Eigen::Vector3d::Ones()*cfg->lookup("acc_error_std")* time_diff;
                stds.segment<3>(0)+=Eigen::Vector3d::Ones()*cfg->lookup("ar_error_std")* time_diff;
            }

            //Convert to stiffness matrix
            stds=stds.cwiseInverse();
        
               //Dynamic cost function 
             ceres::CostFunction *cost_function =
                new ceres::AutoDiffCostFunction<BoulderDynamic, DOF, STATE_SIZE, STATE_SIZE, 6>{
                    new BoulderDynamic{imu->getAngularVelocity(), imu->getAcceleration(), time_diff,stds}};
            problem.AddResidualBlock(cost_function, NULL, &x[STATE_SIZE * t],
                                     &x[STATE_SIZE * (t + 1)], bias);
            
            //Apply ZUPT when grip is true
            if (is_grip[t])
            {
                cost_function =
                    new ceres::AutoDiffCostFunction<ZeroVelocityCost, 3, STATE_SIZE>{new ZeroVelocityCost{}};
                problem.AddResidualBlock(cost_function, NULL, &x[STATE_SIZE * t]);
                (getVel(&x[STATE_SIZE * t], double, )) = VecD<3>::Zero();
            }
            //calculate initial guess via dynamic model ...
            dynamicModel<double>(&x[STATE_SIZE * t], bias, imu->getAngularVelocity(), imu->getAcceleration(), time_diff).toPointer(&x[STATE_SIZE * (t + 1)]);
            double outputQuat[4],outputError[2];
            // ... and orientation estimator
            orient_guess.updateCsgOriEstIMU(imu->getAcceleration().data(),imu->getAngularVelocity().data(),imu->getMagnetometer().data(),outputQuat,magnetometer_bias,outputError);
            adekf::SO3d{outputQuat[0],outputQuat[1],outputQuat[2],outputQuat[3]}.toPointer(&x[STATE_SIZE * (t + 1)]); //storage order of eigen and CCsg is different. This constructor of SO3d converts. 
            
            t++;
            // add magnet cost
            cost_function =
                new ceres::AutoDiffCostFunction<MagnetometerCost, 3, STATE_SIZE,3>{new MagnetometerCost{imu->getMagnetometer()}};
            problem.AddResidualBlock(cost_function, NULL, &x[STATE_SIZE * t],magnetometer_bias);
            
            //plot some raw measurements if needed
            //adekf::viz::plotVector(imu->getMagnetometer(),"Magnet measurement", rowCount,"xyz");
            //adekf::viz::plotVector(imu->getAngularVelocity(),"Angular Velocity", rowCount,"xyz");
            //adekf::viz::plotVector(imu->getAcceleration(),"Acceleration", rowCount,"xyz");
        };
        //add callback to IMU plugin
          imu->addDataCallback<zavi::plugin::FileIMUPlugin, void>(NULL, add_costs);
          //Add bias cost
        ceres::CostFunction *cost_function =
                new ceres::AutoDiffCostFunction<BiasCost, 3, 3>{
                    new BiasCost{}};
            problem.AddResidualBlock(cost_function, NULL,  magnetometer_bias);
    
        LOG(INFO) << "Adding cost functions" << std::endl;

        //run all through data of the imu to build all cost functions, run it fast to avoid long waits
        imu->run(5000, []()
                 { return true; });

        //LOG(INFO) << "Start state: " << getOrient(&x[STATE_SIZE * 1]) LOG_END;
        LOG(INFO) << "IMU Count: " << rowCount LOG_END;
        LOG(INFO) << "Used states: " << t LOG_END;
        LOG(INFO) << "Actual x size: " <<  ((int)rowCount) +7 LOG_END;

        //Add local parametrization to avoid destruction of the orienation quaternion
        ceres::LocalParameterization *localParameterization = STATE_TYPE::getLocalParameterization();
        for (size_t i = 0; i <= t; i++)
        {
            problem.AddParameterBlock(&x[STATE_SIZE * i], STATE_SIZE, localParameterization);
        }

        // Run the solver!
        ceres::Solver::Summary summary;
        LOG(INFO) << "Start_orient_CCSG: " << adekf::SO3d{&x[STATE_SIZE]}  LOG_END;
        LOG(INFO) << "Starting Solver" LOG_END;
        ceres::Solve(options, &problem, &summary);
        std::cout <<  summary.BriefReport() << "\n";
        

        LOG(INFO) << "Start_orient_Batch: " << adekf::SO3d{&x[STATE_SIZE]}  LOG_END;   
        LOG(INFO) << "Magnetometer Bias: " << Eigen::Vector3d(magnetometer_bias).transpose() LOG_END;
        //Calculate covariance
        ceres::Covariance::Options cov_options;
        cov_options.sparse_linear_algebra_library_type=ceres::SUITE_SPARSE;
        ceres::Covariance covariance(cov_options);
        std::vector<std::pair<const double*, const double*> > covariance_blocks;
        //Add all wanted blocks (only covariances for each state x_k and not inter state)
         for (size_t i = 0; i <= t; i++)
        {
            covariance_blocks.push_back(std::make_pair(&x[STATE_SIZE * i], &x[STATE_SIZE * i]));
        }
        bool cov_computation_success=covariance.Compute(covariance_blocks,&problem);
        adekf::aligned_vector<std::tuple<VecD<3>,Eigen::Matrix<double,6,6>,adekf::SO3d> > states_and_covariances;
        //if covariance computation fails default to standard cov 
        if (!cov_computation_success)
            LOG(INFO) << "Defaulting to standard velocity cov" LOG_END;
        //base covariance since the outputted covariance is too low for short transitions
        Eigen::Matrix<double,6,6> standard_cov=Eigen::Matrix<double,6,6>::Identity() * pow((double) cfg->lookup("standard_velocity_std"),2);
        standard_cov.block<3,3>(0,0)=Eigen::Matrix3d::Identity()*pow((double) cfg->lookup("standard_orient_std"),2);
        for (size_t i = 0; i <= t; i++)
        {
            //plot corrected magnetometer
            adekf::SO3d imu_in_world=getOrient(&x[STATE_SIZE * i]);
            adekf::viz::plotVector((imu_in_world*accelerations[i].normalized()-MagnetometerCost::mag_direction).eval(),"Corrected Magnetometer", rowCount,"xyz");
      
            // plot velocity and grip
            VecD<4> to_plot;
            Eigen::Vector3d velocity = getVel(&x[STATE_SIZE * i], double, );
            
            to_plot << velocity, is_grip[i];
            adekf::viz::plotVector(to_plot, "estimated velocity and grip", t + 1, "xyzg");
            //Read covariance blocks
            if(cov_computation_success){
            double cov_pointer[DOF*DOF];
            covariance.GetCovarianceBlockInTangentSpace(&x[STATE_SIZE * i],&x[STATE_SIZE * i],cov_pointer);
            Eigen::Matrix<double,DOF,DOF> state_covariance{cov_pointer};
            states_and_covariances.push_back(std::make_tuple(velocity,standard_cov+state_covariance,imu_in_world));
            }
            else{
                states_and_covariances.push_back(std::make_tuple(velocity,standard_cov,imu_in_world));
                
            }
        }
        LOG(INFO) << "Covariances computed" LOG_END;
        delete[] x;
        //reset IMU to be reusable
        imu->callbacks.clear();
        imu->setTimeRange(start_time,end_time);
        LOG(INFO) << "Finished Batch Estimation" LOG_END;
        return states_and_covariances;
    }
} // velocity batch

#endif // ZAVI_BOULDERN_BATCHESTIMATOR_H
