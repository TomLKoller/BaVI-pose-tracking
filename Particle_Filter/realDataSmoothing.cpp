//
// Created by tomlucas on 03.11.20.

// Enable or disable visualization
#define SHOW_VIZ
#include "ManifoldTypes.hpp"
#include "adekf_viz.h"
#include "ADEKFUtils.h"
#include "TransformIt.h"
#include <filesystem>
#include <nlohmann/json.hpp>
#include "CSV_Reader.hpp"
#include "ZaviConfig.hpp"
#include "datasetFormat.hpp"
#include "ZaVI_Utils.hpp"
#include "Plugins/FileIMUPlugin.hpp"
#include "Plugins/ARTPlugin.hpp"
#include "Plugins/SynchronizedSensorsPlugin.hpp"
#include "Plugins/state_sensor_plugin_posereader.hpp"
#include "Plugins/ResultLoader.hpp"
#include "ClusterParticleFilter.h"
#include "BatchEstimator.h"
#include "probabilistics/mahalonobis_probability.hpp"
#include "gripDetector.h"

namespace fs = std::filesystem;
using json = nlohmann::json;
using namespace adekf::viz;
using namespace Eigen;

// State of the Particle filter
ADEKF_MANIFOLD(Pose3DYaw, , (3, position), (1, yaw))

template <typename T>
using State = Pose3DYaw<T>;
using Cov = Pose3DYawCov;

/**
 * @brief Functor for dynamic model of the PF
 *
 */
struct boulder_dynamic_velocity_only
{
    /**
     * @brief Dynamic model with known world velocity and yaw error sampling
     *
     * @tparam T Scalar Type
     * @tparam NOISE_TYPE noise vector type
     * @param state The state of the particle
     * @param noise A noise pertubation vector with size 7
     * @param input a tuple that contains the world velocity, orientation and time diff
     */
    template <typename T, typename NOISE_TYPE>
    void operator()(State<T> &state, const NOISE_TYPE &noise, const std::tuple<Eigen::Vector3d, adekf::SO3d, double> input) const
    {
        double deltaT = std::get<2>(input);
        // sample the yaw
        state.yaw += NOISE(6, 1) * deltaT;
        Eigen::Matrix<T, 3, 1> yaw_vector;
        yaw_vector << T(0.), T(0.), state.yaw;
        // Yaw error quaternion
        adekf::SO3<T> yaw_offset{yaw_vector};
        // Update position by rotated velocity
        state.position += (yaw_offset * (std::get<0>(input) + NOISE(3, 3))) * deltaT;
    }
} boulder_dynamic_velocity_only_;

/**
 * @brief Reads the handles positions from files
 *
 * @param handle_file file path of the handles
 * @return adekf::aligned_vector<Eigen::Vector3d> A vector with all handle center positions
 */
adekf::aligned_vector<Eigen::Vector3d> readHandles(std::string handle_file)
{
    // load file
    std::ifstream file{handle_file};
    adekf::aligned_vector<Eigen::Vector3d> handles;
    if (file.good())
    {
        std::string dummy;
        // skip header line
        getline(file, dummy);
    }
    // read all positions from csv
    while (file.good())
    {
        Eigen::Vector4d handle;
        if (zavi::csv_reader::readLineToEigen(file, handle, ','))
        {
            // only read position not index
            handles.push_back(handle.segment<3>(1));
        }
    }

    return handles;
}

/**
 * @brief Read handles of a route from a folder
 *
 * @param map_folder The folder to search the data in
 * @param trackID The route ID
 * @return adekf::aligned_vector<Eigen::Vector3d> A vector with all handles of the route
 */
adekf::aligned_vector<Eigen::Vector3d> readHandlesFromFolder(const std::string map_folder, const std::string &trackID)
{
    LOG(INFO) << "Handle File: " << map_folder + "/routes/" + trackID + ".csv" LOG_END;
    return readHandles(map_folder + "/routes/" + trackID + ".csv");
}

/**
 * @brief Computes the most likely path of a Hand  IMU given the route handle positions.
 *
 * Implements the algorithm described in our Paper
 *
 * Event-Domain Knowledge in Inertial Sensor Based State Estimation of Human Motion
 * Tom L. Koller, Tim Laue and Udo Frese
 *
 *
 * @param imu The IMUPlugin with the hand data
 * @param art The ArtPlugin with the GT data (for visualization purpose only)
 * @param handles A vector with the center positions of the route holds
 * @param view_handles A vector with the center positions of all holds on the wall
 * @param start_handles The possible start holds of the route
 * @param output_path The output path for the Particle Filter
 * @param smooth_output_path The output path for the Particle Smoother
 * @param transition_output_path The output path for the Transition Estimator
 * @param trial_data The trial_data.json data
 * @param magnet_in_art The magnet direction in world coordinates
 * @param known_start Use start_handles or localize from all holds of the route
 */
void calcHand(std::shared_ptr<zavi::plugin::FileIMUPlugin> imu, std::shared_ptr<zavi::plugin::ARTPlugin> art, const adekf::aligned_vector<Eigen::Vector3d> &handles, const adekf::aligned_vector<Eigen::Vector3d> &view_handles, const size_t start_handles[2], const std::string &output_path, const std::string &smooth_output_path, const std::string &transition_output_path, const bavi::json &trial_data, const Eigen::Vector3d &magnet_in_art, bool known_start = true)
{
    // Reset viz
    adekf::viz::PoseRenderer::removeAllActors();
    State<double> start_state;
    // read times from trial json
    double start_time = trial_data["start_time"];
    double descent_time = trial_data["descent_time"];
    art->setTimeRange(start_time);

    Cov start_cov = start_cov.Identity() * 0.001;
    // cov of state.yaw
    start_cov(3, 3) = pow(cfg->lookup("yaw_start_std"), 2);
    // Collect possible start positions from route
    adekf::aligned_vector<State<double>> start_states;
    start_state.yaw = Eigen::Matrix<double, 1, 1>{0.};
    if (known_start)
    {
        for (size_t i = 0; i < 2; i++)
        {
            start_state.position = handles[start_handles[i]];
            start_states.push_back(start_state);
        }
    }
    else
    {
        for (auto handle : handles)
        {
            start_state.position = handle;
            start_states.push_back(start_state);
        }
    }
    size_t num_particles = cfg->lookup("num_particles");
    // Adjust number of particles to be evenly distributed among handles. Adds at most n-1 particles.
    num_particles += num_particles % start_states.size();
    // INitialize the PF
    auto pf = std::make_shared<ClusterParticleFilter<State<double>>>(num_particles, start_states, start_cov);
    // Stride for particle visualization (only visualizes every strideth particle)
    int stride = cfg->lookup("particle_stride");
    PoseRenderer::displayPoints(handles, "blue", 0.03);
    // Particle viz object
    auto pc = PoseRenderer::displayParticles(transformToPoints(
                                                 pf->particles_, [](auto state)
                                                 { return state.position; },
                                                 stride),
                                             "green", 0.02);

    Eigen::Matrix3d grip_cov = grip_cov.Identity() * pow((double)cfg->lookup("grip_std"), 2);
    /**
     * @brief Probability function for the position at gripping
     * @param particle the particle state
     * @param handles the hold center positions
     *
     */
    auto grip_set_probability_field = [&grip_cov](auto particle, auto handles)
    {
        // struct that can compute the probability of a Set-Gaussian mixture  with uniform minimum
        static mahalonobis_probability_with_set measurement_probability{grip_cov, cfg->lookup("grip_set"), cfg->lookup("max_std_for_set")};
        ScalarOf(particle) prob_sum{0.};
        // Compute the highest likelihood for all gripping positions
        for (auto handle : handles)
        {
            prob_sum = std::max(prob_sum, exp(measurement_probability(particle.position, handle)));
        }
        return prob_sum;
    };

    /**
     * @brief Function to cluster the data based on probability to belong to a position
     *
     */
    auto get_best_cluster = [&grip_cov](auto particle, auto handles)
    {
        // Object to calculate probability
        static mahalonobis_probability_with_set measurement_probability{grip_cov, cfg->lookup("grip_set"), cfg->lookup("max_std_for_set")};
        // initialize best position to elsewhere
        ScalarOf(particle) max_prob{measurement_probability.min_return};
        size_t best_cluster = 0; // Else class
        for (size_t index = 0; index < handles.size(); index++)
        {
            double prob = exp(measurement_probability(particle.position, handles[index]));
            if (prob > max_prob)
            {
                best_cluster = index + 1;
                max_prob = prob;
            }
        }
        return best_cluster;
    };

    LOG(INFO) << "Trial starts at: " << start_time LOG_END;

    // Find all grips
    imu->setTimeRange(start_time, trial_data["end_time"]);
    // The start index of the time window in IMU indices
    size_t imu_current_row = imu->getCurrentRow();
    auto is_grip = gripDetector(imu, cfg->lookup("std_thresh"));
    LOG(INFO) << "Grip vector size: " << is_grip.size() LOG_END;

    // Find last grip
    size_t last_grip = 0;
    size_t first_grip = is_grip.size();
    for (size_t index = 0; index < is_grip.size(); index++)
    {
        if (is_grip[index])
        {
            last_grip = index;
            if (index < first_grip)
            {
                first_grip = index;
            }
        }
    }
    size_t t = 0;
    // last grip starts from the beginning of the time window.
    double end_time = imu->getTimeAt(imu_current_row + last_grip);

    LOG(INFO) << "Start time: " << start_time << " Descent Time: " << descent_time << " End Time: " << end_time LOG_END;
    // Reset imu to shortened window
    imu->setTimeRange(start_time, end_time);

    // Create output file
    std::ofstream output{output_path};
    output << "Time in s, x in m, y in m, z in m, imu_in_world w,imu_in_world x,imu_in_world y,imu_in_world z,vx in m/s,vy in m/s,vz in m/s,cluster id,cluster probability" << std::endl;

    // run Transition Estimation
    auto state_and_covariance = velocity_batch::transitionEstimation(imu, start_time, end_time, is_grip, magnet_in_art);
    // Setup dynamic cov (static part)
    Eigen::Matrix<double, 7, 7> dynamic_cov = dynamic_cov.Zero();
    dynamic_cov(6, 6) = pow(cfg->lookup("yaw_std"), 2);
    // Vector to store inputs, covariances and time stamps (for smoothing)
    adekf::aligned_vector<std::tuple<Eigen::Vector3d, adekf::SO3d, double>> inputs;
    adekf::aligned_vector<Cov> covs;
    std::vector<double> times;
    /**
     * @brief IMU callback that runs the PF on every IMU sample
     *
     */
    auto imu_callback = [&get_best_cluster, &times, &covs, &inputs, first_grip, &dynamic_cov, &start_time, &view_handles, &descent_time, &handles, &output, &pc, &is_grip, &t, state_and_covariance, &grip_set_probability_field, &stride](decltype(imu.get()) imu, decltype(pf.get()) pf, double deltaT)
    {
        // Store time stamps and inputs
        times.push_back(imu->nextTime);
        inputs.emplace_back(std::get<0>(state_and_covariance[t]), std::get<2>(state_and_covariance[t]), deltaT);
        // Do not start PF before first grip
        if (t < first_grip)
        {
            t++;
            return;
        }
        LOG(INFO) << "IMU callback at: " << imu->nextTime LOG_END if (t >= state_and_covariance.size() || t >= is_grip.size())
        {
            LOG(WARNING) << "Trying to read unavailable velocity" LOG_END;
            return;
        }
        // read cvoariance
        dynamic_cov.block<6, 6>(0, 0) = std::get<1>(state_and_covariance[t]);
        // Compute the dynamic covariance for smoothing (transform non-additive to additive)
        auto mean = pf->getMean();
        constexpr size_t DOF = adekf::DOFOf<State<double>>;
        constexpr size_t NOISE_DIM = 7;
        // auto derivation of the dynamic model
        auto derivator = adekf::getDerivator<DOF + NOISE_DIM>();
        auto input = eval(mean + derivator.template head<DOF>());
        boulder_dynamic_velocity_only_(input, derivator.template tail<NOISE_DIM>(), inputs.back());
        // Jacobian of the dynamic model
        auto F = adekf::extractJacobi(input - mean);
        covs.push_back(F.template rightCols<NOISE_DIM>() * dynamic_cov * F.template rightCols<NOISE_DIM>().transpose());

        pf->predict(boulder_dynamic_velocity_only_, dynamic_cov, inputs.back());

        // Perform grip update
        if (is_grip[t])
        {
            // Use route map during ascent
            if (descent_time < start_time || imu->nextTime < descent_time)
            {
                pf->update(grip_set_probability_field, handles);
                pf->updateClusters(handles, get_best_cluster);
            }
            else
            {
                // use Wall map during descent
                pf->update(grip_set_probability_field, view_handles);
                pf->updateClusters(view_handles, get_best_cluster);
            }
        }
        LOG(INFO) << "Starting resampling" LOG_END;
        pf->resample(true);
        LOG(INFO) << "Storing PDF" LOG_END;
        pf->storePDF();

        LOG(INFO) << "Writing Output" LOG_END

            /**
             * @brief Function to print a cluster mean to a string
             * @param mean the particle mean
             * @param id the cluster id
             * @param weight the cluster weight
             */
            auto output_mean = [&state_and_covariance, &imu, &t](auto mean, int id, double weight)
        {
            auto position = mean.position;
            Eigen::Matrix<double, 3, 1> yaw_vector;
            yaw_vector << 0., 0., mean.yaw;
            adekf::SO3<double> yaw_offset{yaw_vector};

            auto velocity = yaw_offset * std::get<0>(state_and_covariance[t]);
            adekf::SO3<double> orient = yaw_offset * std::get<2>(state_and_covariance[t]);
            LOG(INFO) << mean;
            std::stringstream out;
            out << imu->nextTime << "," << position.x() << "," << position.y() << "," << position.z() << "," << orient.w() << "," << orient.x() << "," << orient.y() << "," << orient.z() << "," << velocity.x() << "," << velocity.y() << "," << velocity.z() << "," << id << "," << weight << std::endl;
            ;
            return out.str();
        };
        mean = pf->getMean();
        output << output_mean(mean, -1, 1);
        LOG(INFO) << "Writing Clusters" LOG_END;
        std::map<size_t, std::vector<size_t>> clusters;
        // Collect particles of all clusters
        for (size_t i = 0; i < pf->clusters_.size(); i++)
        {
            size_t nearest = pf->clusters_[i];
            if (clusters.find(nearest) == clusters.end())
            {
                clusters.emplace(nearest, std::vector<size_t>{});
            }
            clusters[nearest].push_back(i);
        }
        // Output cluster means
        for (auto pair : clusters)
        {
            mean = pf->getMean(pair.second.begin(), pair.second.end(), pf->particles_, pf->weights_);
            output << output_mean(mean, pair.first, pf->getWeightSum(pair.second.begin(), pair.second.end(), pf->weights_));
        }
// Update visualization
#ifdef SHOW_VIZ
        auto points = transformToPoints(
            pf->particles_, [](auto state)
            { return state.position; },
            stride);
        for (int i = 0; i < pc->GetNumberOfPoints(); i++)
        {
            pc->InsertPoint(i, points[i].data());
        }
        pc->Modified();
#endif
        t++;
    };

    // Read start for transition integration
    State<double> transition_integration;
    transition_integration.yaw << 0.;
    bool has_started = false;
    size_t start_t;
    /**
     * @brief Callback to read first art measurement for transition integration
     *
     */
    auto art_callback = [&has_started, &t, &start_t, &transition_integration](decltype(art.get()) art, decltype(pf.get()) pf, double deltaT)
    {
        if (!has_started)
        {
            has_started = true;
            start_t = t;
            transition_integration.position = art->getPosition();
        }
    };
    // Add callbacks
    art->addDataCallback(art, pf, art_callback);
    imu->addDataCallback(imu, pf, imu_callback);
    // synchronize sensor reading
    auto sync = zavi::plugin::SynchronizedSensorPlugin({art, imu}, false);

    // Viz object of the ground truth
    auto gt_render = adekf::viz::PoseRenderer::displayPoseGeneric<zavi::plugin::ARTPlugin, zavi::plugin::SensorPoseDelayedReader>(art.get(), "red", 0.1);
    // Functor to check whether the IMU should go on running
    auto RUNNING = []()
    { return !adekf::viz::PoseRenderer::isDone(); };

    // Run sensors with configured frequency
    sync.run(cfg->lookup("loop_frequency"), RUNNING);
    LOG(INFO) << "Used up: " << t << "values of is_grip" LOG_END;
    /**
     * @brief Functor to output smoothed estimate
     *
     */
    auto output_mean = [&state_and_covariance](auto mean, int id, double weight, double time, size_t orient_index)
    {
        auto position = mean.position;
        Eigen::Matrix<double, 3, 1> yaw_vector;
        yaw_vector << 0., 0., mean.yaw;
        adekf::SO3<double> yaw_offset{yaw_vector};
        auto velocity = yaw_offset * std::get<0>(state_and_covariance[orient_index]);

        adekf::SO3<double> orient = yaw_offset * std::get<2>(state_and_covariance[orient_index]);
        std::stringstream out;
        out << time << "," << position.x() << "," << position.y() << "," << position.z() << "," << orient.w() << "," << orient.x() << "," << orient.y() << "," << orient.z() << "," << velocity.x() << "," << velocity.y() << "," << velocity.z() << "," << id << "," << weight << std::endl;
        ;
        return out.str();
    };

    // Run smoother
    pf->smooth(boulder_dynamic_velocity_only_, covs, inputs, first_grip);
    // OPen smooth output file
    std::ofstream smooth_output{smooth_output_path};
    smooth_output << "Time in s, x in m, y in m, z in m, imu_in_world w,imu_in_world x,imu_in_world y,imu_in_world z,vx in m/s,vy in m/s,vz in m/s,cluster id,cluster probability" << std::endl;

    // Output smoothed data
    for (size_t old_index = 0; old_index < pf->old_pdfs_.size(); old_index++)
    {
        auto pdf = pf->old_pdfs_[old_index];
        // Output mean
        auto range = boost::irange(size_t(0), pdf.weights_.size());
        auto mean = pf->getMean(range.begin(), range.end(), pdf.particles_, pdf.weights_);
        smooth_output << output_mean(mean, -1, 1, times[first_grip + old_index], first_grip + old_index);

        // Find clusters
        std::map<size_t, std::vector<size_t>> clusters;
        for (size_t i = 0; i < pdf.clusters_.size(); i++)
        {
            size_t nearest = pdf.clusters_[i];
            if (clusters.find(nearest) == clusters.end())
            {
                clusters.emplace(nearest, std::vector<size_t>{});
            }
            clusters[nearest].push_back(i);
        }
        // Output cluster means
        for (auto pair : clusters)
        {
            auto mean = pf->getMean(pair.second.begin(), pair.second.end(), pdf.particles_, pdf.weights_);
            smooth_output << output_mean(mean, pair.first, pf->getWeightSum(pair.second.begin(), pair.second.end(), pdf.weights_), times[first_grip + old_index], first_grip + old_index);
        }
    }
    /**
     * @brief Functor to output transition estimation
     *
     */
    auto output_transition = [&state_and_covariance](auto mean, double time, size_t orient_index)
    {
        auto position = mean.position;
        Eigen::Matrix<double, 3, 1> yaw_vector;
        yaw_vector << 0., 0., mean.yaw;
        adekf::SO3<double> yaw_offset{yaw_vector};

        auto velocity = yaw_offset * std::get<0>(state_and_covariance[orient_index]);
        adekf::SO3<double> orient = yaw_offset * std::get<2>(state_and_covariance[orient_index]);
        std::stringstream out;
        out << time << "," << position.x() << "," << position.y() << "," << position.z() << "," << orient.w() << "," << orient.x() << "," << orient.y() << "," << orient.z() << "," << velocity.x() << "," << velocity.y() << "," << velocity.z() << std::endl;
        ;
        return out.str();
    };
    // Open transition output
    std::ofstream transition_output{transition_output_path};
    transition_output << "Time in s, x in m, y in m, z in m, imu_in_world w,imu_in_world x,imu_in_world y,imu_in_world z,vx in m/s,vy in m/s,vz in m/s" << std::endl;
    // Write start value
    transition_output << output_transition(transition_integration, times[start_t] - 0.01, start_t);
    // Output result of the dynamic model to transition
    for (size_t i = start_t; i < t; i++)
    {
        if (i > times.size())
        {
            LOG(WARNING) << "Transition out of bounds" LOG_END;
        }
        boulder_dynamic_velocity_only_(transition_integration, Eigen::Matrix<double, 7, 1>::Zero(), inputs[i]);
        transition_output << output_transition(transition_integration, times[i], i);
    }

    // Close all files
    output.flush();
    output.close();
    smooth_output.flush();
    smooth_output.close();
    transition_output.flush();
    transition_output.close();

    LOG(INFO) << "Finished one PF" LOG_END;
    adekf::viz::PoseRenderer::removeActor(gt_render);
}

/**
 * @brief Runs the PF Smoother on the whole dataset
 *
 * @param argc
 * @param argv
 * @return int
 */
int main(int argc, char **argv)
{
    if (argc != 7)
    {
        std::cout << "Passed: ";
        for (int i = 0; i < argc; i++)
            std::cout << argv[i] << " ";
        std::cout LOG_END;
        std::cout << "Expected Format: Zavi-Particle_Filter_Smooth <dataset_directory>  <result_directory> <smooth_result_directory> <transition_result_directory> <map_directory> <Config_folder>" << std::endl;
        return -1;
    }
    adekf::viz::initGuis(argc, argv);
    // reset locale for csv reader
    setlocale(LC_ALL, "C");
    printf("Starting up\n");
    initConfig(argv[6], "main_config.cfg");
    // Setup logging
    google::InitGoogleLogging(argv[0]);
    google::SetStderrLogging(cfg->lookup("log_level"));
    // init timer
    zavi::Timer::initTime(cfg->lookup("time_lapse"));
    // Read filter values in config to only evaluate a subset of the dataset
    std::string target_trial = cfg->lookup("target_trial");
    std::string target_sensor = cfg->lookup("target_sensor");
    if (target_trial != "_")
    {
        LOG(INFO) << "Only going for Trial: " << target_trial LOG_END;
    }
    // load wall maps
    adekf::aligned_vector<Eigen::Vector3d> v01_handles, v02_handles;
    std::string v01_routes[] = {"yellow1", "blue1", "black1", "green1", "purple1", "frog1", "shell1"};
    for (auto name : v01_routes)
    {
        auto handles = readHandlesFromFolder(argv[5], name);
        v01_handles.insert(v01_handles.end(), handles.begin(), handles.end());
    }
    std::string v02_routes[] = {"black2", "blue2", "frog2", "purple2", "red2", "white2", "yellow2"};
    for (auto name : v02_routes)
    {
        auto handles = readHandlesFromFolder(argv[5], name);
        v02_handles.insert(v02_handles.end(), handles.begin(), handles.end());
    }
    // Open second thread for PF smoother to enable visualization
    std::thread loop{
        [&]()
        {
                        std::filesystem::path dataset_path{argv[1]};
                        //Iterate participants
                        for (const auto &participant_entry : bavi::iterateParticipants(dataset_path))
                        {
                         std::filesystem::path	participant_folder{participant_entry.path()};
	                     participant_folder = participant_folder.lexically_normal();
                         std::cout << participant_folder <<std::endl;
    
                         std::filesystem::path participant_data{participant_folder.string()+ "/participant_data.json"};
                         participant_data = participant_data.lexically_normal();
                         //Copy participant_data.json to output paths
                        for(char * folder: {argv[2],argv[3],argv[4]}){
                         std::filesystem::path participant_data_output = std::string(folder) + (participant_data).lexically_relative(dataset_path).string();
                         fs::create_directories(participant_data_output.parent_path());
                         std::filesystem::copy_file(participant_data, participant_data_output, std::filesystem::copy_options::update_existing);
                        }
                        //Read magnet direction
                        auto dataset_json=bavi::readJson(dataset_path,"data_set.json");
                        auto magnet_in_art_array1 = dataset_json["magnet_in_art_V01"];
                        Eigen::Vector3d magnet_in_art1{magnet_in_art_array1[0], magnet_in_art_array1[1], magnet_in_art_array1[2]};
                        auto magnet_in_art_array2 = dataset_json["magnet_in_art_V02"];
                        Eigen::Vector3d magnet_in_art2{magnet_in_art_array2[0], magnet_in_art_array2[1], magnet_in_art_array2[2]};
                            

                         // Iterate all viewpoints
                         for (const auto &view_folder : bavi::iterateViews(participant_folder))
                         {
                             //Select wall map
                             bool is_v01 = view_folder.path().string().find("V01") != std::string::npos;
                             Eigen::Vector3d view_magnet;
                             if (is_v01)
                                view_magnet=magnet_in_art1;
                            else
                                view_magnet=magnet_in_art2;
                             // Iterate all trials
                             for (const auto &trial_folder : bavi::iterateTrials(view_folder.path()))
                             {  
                                //Skip if not the target trial
                                 if (target_trial != "_" && trial_folder.path().string().find(target_trial)==std::string::npos)
				                        continue; 
                                std::cout << trial_folder LOG_END;
    
                                 // Copy trial_data.json to output folders
                                std::filesystem::path trial_data_file = trial_folder.path().string() + "/trial_data.json";
                                 trial_data_file = trial_data_file.lexically_normal();
                                 for(char * folder: {argv[2],argv[3],argv[4]}){
                        
                                 std::filesystem::path trial_data_output = std::string(folder) + std::filesystem::path(trial_data_file).lexically_relative(dataset_path).string();
                                 fs::create_directories(trial_data_output.parent_path());
                                 std::filesystem::copy_file(trial_data_file, trial_data_output, std::filesystem::copy_options::overwrite_existing);
                                 }
                                 //Load trial data
                                 auto trial_data = bavi::readJson(trial_folder, "trial_data.json");
                                 LOG(INFO) << "Loading handles." LOG_END;
                                 auto handles = readHandlesFromFolder(argv[5], trial_data["TrackID"]);
                                 std::filesystem::path map_folder{argv[5]};
                                 //Read start handles
                                 auto start_handles_json = bavi::readJson(map_folder, "start_handles.json");
                                 size_t start_handles[2];
                                 nlohmann::json start_handle_array = start_handles_json[(std::string)trial_data["TrackID"]];
                                 start_handles[0] = start_handle_array[0];
                                 start_handles[1] = start_handle_array[1];
                                
                                 // Skip if no video available since real trial start time is unknown in this cases
                                 if (trial_data["video_synched"] == 0)
                                 {
                                     LOG(INFO) << "Skipped due to missing video (no start_time available)" LOG_END; continue;
                                 }

                                 // Find all hand imu files
                                 for (const auto &imu_file : bavi::iterateFolder(trial_folder, std::regex{"imu_\\w+hand.csv"}, bavi::matchesRegex))
                                 {
                                     //Skip if not the target sensor
                                     if (target_sensor != "_" && imu_file.path().string().find(target_sensor)==std::string::npos)
					                    continue; 
                                    //find sensor name
                                    std::regex fex{"(.+?)(imu_)(\\w+)(.csv)"};
                                    std::smatch match;
                                    std::string temp=imu_file.path().string();
                                    //std::cout << temp LOG_END;
                                     if(std::regex_search(temp,match, fex)){
                                         std::string sensor_name =match[3];
                                         //std::cout<< sensor_name LOG_END;
                                         // skip if no data is avaialbe for this sensor (for example broken record  or sensor rip off)
                                        if (trial_data["data_mask"][sensor_name]==0) 
                                            continue;
                                     }
                                     //Load IMU data
                                     auto imu = std::make_shared<zavi::plugin::FileIMUPlugin>(imu_file.path().string(), 0, -1, false);
                                     imu->loadAll();
                                     //Load ground truth data (for visualization)
                                     std::regex ex{"(imu_)(\\w+.csv)"};
                                     std::string art_name = std::regex_replace(imu_file.path().string(), ex, "art_world_$2");
                                     
                                    
                                    auto art = std::make_shared<zavi::plugin::ARTPlugin>(art_name.c_str(), std::set<int>{});
                                     LOG(INFO) << "Art Path:" << art_name << std::endl;
                                     //Create output paths
                                     std::string output_name = std::regex_replace(imu_file.path().string(), ex, "pose_$2");
                                     std::string output_path = std::string(argv[2]) + std::filesystem::path(output_name).lexically_relative(dataset_path).string();
                                     std::string smooth_output_path=   std::string(argv[3]) + std::filesystem::path(output_name).lexically_relative(dataset_path).string();
                                     std::string transition_output_path=   std::string(argv[4]) + std::filesystem::path(output_name).lexically_relative(dataset_path).string();
                                     
                                     LOG(INFO) << "Output path: " << output_path LOG_END;
                                     //run position estimator
                                     calcHand(imu, art, handles, is_v01 ? v01_handles : v02_handles, start_handles, output_path, smooth_output_path,transition_output_path,  trial_data, view_magnet, cfg->lookup("known_start_handles"));
                                 }
                             }
                         }
                     } }};

#ifdef SHOW_VIZ
    LOG(INFO) << "Starting Gui" LOG_END
        adekf::viz::runGuis();
#endif
    zavi::ThreadedObject::joinAll();
    loop.join();
    return 0;
}