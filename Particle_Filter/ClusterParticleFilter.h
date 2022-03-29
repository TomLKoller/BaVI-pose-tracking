//
// Created by Tom Koller on 03.11.20.
//

#ifndef ZAVI_BOULDERN_CLUSTERPARTICLEFILTER_H
#define ZAVI_BOULDERN_CLUSTERPARTICLEFILTER_H

#include <vector>
#include <Eigen/Core>
#include "MultiVariateGaussianNoiseVector.h"
#include "ADEKFUtils.h"
#include "ZaviConfig.hpp"
#include <algorithm>
#include <boost/range/irange.hpp>


/**
 * @brief Particle Filter with clustering capability
 * Based on: A. Doucet and A. M. Johansen, “A tutorial on particle filtering and
 * smoothing: Fifteen years later,” Handbook of nonlinear filtering, no. 12,
 * pp. 656–704, 2009.   
 * 
 * @tparam State_ The state description (Eigen::Matrix or adekf::Manifold)
 */
template <class State_>
class ClusterParticleFilter
{
private:
    /**
     * @brief Struct to store a pdf
     * 
     */
    struct PDF
    {
        adekf::aligned_vector<State_> particles_;
        std::vector<double> weights_;
        std::vector<size_t> clusters_;
    };

public:
    //The state type
    using State = State_;
    //All particles' states
    adekf::aligned_vector<State> particles_;
    //Particle indexes
     std::vector<size_t> clusters_;
     //Particle weights
    std::vector<double> weights_;
    //Old pdfs for smoothing
    std::vector<PDF> old_pdfs_;
    //random engine for prediction
    static inline std::default_random_engine generator_;
    //normal distribution to sample from
    std::uniform_real_distribution<double> distribution_;

    /**
     * @brief Set the Start state of the particle filter as a Gaussian
     * 
     * @param start_state Mean of the Gaussian
     * @param cov  Covariance of the Gaussian
     */
    void setStartState(const State &start_state, const adekf::CovarianceOf<State> &cov){
        adekf::MultiVariateGaussianNoiseVector noise{cov.col(0).Zero(), cov};
        //Add noise to particle
        for (auto &particle : particles_)
        {
            particle = particle + noise.poll();
        }
    }
    /**
     * @brief Construct a new Cluster Particle Filter with a Gaussian start state 
     * 
     * @param n number of particles to use
     * @param start_state Mean of the start state
     * @param cov Covariance of the start state
     */
    ClusterParticleFilter(size_t n, const State &start_state, const adekf::CovarianceOf<State> &cov) : particles_(n, start_state), clusters_(n,0), weights_(n, 1. / n),
                                                                                                distribution_(0, 1. / n)
    {
        setStartState(start_state,cov);
    }

/**
 * @brief Construct a new Cluster Particle Filter object with multi-modal Gaussian start state
 * 
 * @param n number of particles
 * @param start_states Means of the Gaussians
 * @param cov Covariance of all Gaussians
 */
ClusterParticleFilter(size_t n, const adekf::aligned_vector<State> &start_states, const adekf::CovarianceOf<State> &cov) : particles_(n, start_states.front()), clusters_(n,0), weights_(n, 1. / n),
                                                                                                distribution_(0, 1. / n)
    {
        adekf::MultiVariateGaussianNoiseVector noise{cov.col(0).Zero(), cov};
       //Distribute particles uniformly per start state
       size_t particles_per_state=n/start_states.size();
       size_t base_index=0;
       for(size_t j=0; j < start_states.size(); j++){

           for(size_t i=0; i < particles_per_state; i++){
               particles_[base_index+i]=start_states[j]+noise.poll();
               //Assign first clusters
               clusters_[base_index+i]=j;
           }
           base_index+=particles_per_state;
       }
    }

    /**
     * @brief Normalizes sum of weights to 1.
     * 
     * Checks if probability distribution has degenerated (sum of weights =0)
     * 
     * @param weights The weights to normalize
     * @return true If successfull
     * @return false If degenerated
     */
    bool normaliseWeights(std::vector<double> &weights)
    {
        double norm_factor = std::accumulate(weights.begin(), weights.end(), 0.);
        if (norm_factor != 0.)
        {
            norm_factor=1./norm_factor;
            std::transform(weights.begin(), weights.end(), weights.begin(), [&](auto value)
                           { return value * norm_factor; });
                           return true;
        }
        else
        {
            LOG(WARNING) << "Probability distribution degenerated completely. Setting weights uniformly." LOG_END;
            norm_factor=1./weights.size();
            std::transform(weights.begin(), weights.end(), weights.begin(), [&](auto value)
                           { return  norm_factor; });
                           return false;
        }
    }
    /**
     * @brief Shorthand to call normaliseWeights on current weights
     */
    bool normaliseWeights()
    {
        return normaliseWeights(weights_);
    }

    /**
     * @brief Stores the current pdf for smoothing
     * 
     */
    void storePDF()
    {
        old_pdfs_.push_back(PDF{particles_, weights_,clusters_});
    }
    /**
     * @brief Prediction step of the particle Filter
     * 
     * Does not change the weights but only the particle states.
     * 
     * @tparam DYNAMIC_MODEL Class of dynamic Model
     * @tparam COV_TYPE  Type of covariance
     * @tparam ARGS Variadic arguments for the dynamic model
     * @param dynamicModel The dynamic Model that defines the particle transition. Has to change the particle in place. Called as dynamicModel(state,noise vector, args...)
     * @param cov Covariance of the transition. How the noise is applied has to be defined in the dynamic model. 
     * @param args Additional arguments for the dynamic model
     */
    template <class DYNAMIC_MODEL, class COV_TYPE, typename... ARGS>
    void predict(const DYNAMIC_MODEL &dynamicModel, const COV_TYPE &cov, const ARGS &...args)
    {
        adekf::MultiVariateGaussianNoiseVector noise{cov.col(0).Zero(), cov};
        for (size_t i = 0; i < particles_.size(); i++)
        {
            dynamicModel(particles_[i], noise.poll(), args...);
        }
    }
    /**
     * @brief Update step of the particle filter. 
     * 
     * Normalizes weights
     * 
     * @tparam PROBABILITY_MODEL Class of the probabilityModel
     * @tparam ARGS Variadic arguments for the probabilityModel
     * @param probabilityModel the model to compute the particle likelihood. 
     * @param args Additional arguments for the model
     */
    template <typename PROBABILITY_MODEL, typename... ARGS>
    void update(const PROBABILITY_MODEL &probabilityModel, const ARGS &...args)
    {
        for (size_t i = 0; i < particles_.size(); i++)
        {
            weights_[i] *= probabilityModel(particles_[i],args ...);
        }
        normaliseWeights();
    }


    /**
     * @brief Update step of the cluster filter where Gaussian distribution of the measurement is assumed.
     * 
     * Calculates the measurement's Gaussian likelihood to update the particle's weight.
     * Normalizes all weights.
     * 
     * @tparam MEASUREMENT_MODEL Class of the measurement model
     * @tparam MEASUREMENT Class of the measurement
     * @param measurement_model The measurement model h(x)
     * @param measurement The actual measurement
     * @param covariance Covariance of the measurement
     */
    template <typename MEASUREMENT_MODEL, typename MEASUREMENT>
    void updateGaussianMeasurement(const MEASUREMENT_MODEL &measurement_model, const MEASUREMENT &measurement, const adekf::CovarianceOf<MEASUREMENT> &covariance)
    {
        mahalonobis_probability measurement_probability{covariance};
        for (size_t i = 0; i < particles_.size(); i++)
        {
            weights_[i] *= exp(measurement_probability(measurement_model(particles_[i]), measurement));
        }
        normaliseWeights();
    }


    /**
     * @brief Find the closest cluster for a particle based on a distance metric
     * 
     * @tparam CLUSTER_DISTANCE Class of the distance function. 
     * @tparam CLUSTER_MEANS array like list of cluster means
     * @tparam ARGS Variadic arguments for distance function
     * @param particle The particle to find the cluster for
     * @param cluster_means The possible cluster means
     * @param clusterDistance The distance function between cluster mean and particle
     * @param args additional arguments for the cluster function
     * @return size_t The cluster index as it is in cluster_means
     */
    template <typename CLUSTER_DISTANCE,typename CLUSTER_MEANS, typename... ARGS>
    size_t findClosestCluster(const State &particle, const CLUSTER_MEANS & cluster_means, const CLUSTER_DISTANCE &clusterDistance, const ARGS &...args)
{

    double distance = clusterDistance(particle, cluster_means[0], args ...);
    size_t closest = 0;
    for (size_t i = 1; i < cluster_means.size(); i++)
    {
        double t_dist = clusterDistance(particle, cluster_means[i], args ...);
        if (t_dist < distance)
        {
            distance = t_dist;
            closest = i;
        }
    }
    return closest;
};

    /**
     * @brief Updates cluster indices based on a custom clustering function
     * 
     * @tparam BEST_CLUSTER The function to cluster
     * @tparam CLUSTER_MEANS Array like class of cluster means
     * @tparam ARGS Variadic arguments for cluster function
     * @param cluster_means The means of available clusters. Pass empty list if not required.
     * @param get_best_cluster Function that returns the index of the best cluster.
     * @param args Additional arguments for cluster function
     */
    template <typename BEST_CLUSTER,typename CLUSTER_MEANS, typename... ARGS>
    void updateClusters(const CLUSTER_MEANS & cluster_means,const BEST_CLUSTER &get_best_cluster, const ARGS &...args)
    {
        for (size_t i = 0; i < clusters_.size(); i++)
        {
            clusters_[i] = get_best_cluster(particles_[i], cluster_means,  args ...);
        }
    }

    /**
     * @brief Applies a set constraint to the particles. 
     * 
     * Sets particle weight to zero if it does not fulfill the constraint. 
     * Normalizes the weights. 
     * 
     * @tparam CONSTRAIN_MODEL Class of the constraint model
     * @param constrain_model A function that returns true if the particles state is allowed. 
     */
    template <typename CONSTRAIN_MODEL>
    void constrain(const CONSTRAIN_MODEL &constrain_model)
    {
        for (size_t i = 0; i < particles_.size(); i++)
        {
            if (!constrain_model(particles_[i]))
                weights_[i] = 0.;
        }
        normaliseWeights();
    }

    /**
     * @brief Performs stochastic resampling on the particles.
     * 
     * Can use the Effective sample size  criterion (ESS) to decide when to resample. 
     * Sets all weights to 1/n after resampling.
     * 
     * Performs in place resampling.
     * 
     * @param use_ess If true resampling is only performed when the ESS is met. Else resample on every call.
     * @param weights The weights to be resampled. 
     * @param particles The states to be resampled.
     * @param clusters The clusters to be resampled. 
     */
    void resample(bool use_ess, std::vector<double> &weights, adekf::aligned_vector<State> &particles, std::vector<size_t> & clusters)
    {

        if (use_ess)
        {
            // Compute ESS and check if greater than configured threshold
            double ess = 1. / std::accumulate(weights.begin(), weights.end(), 0., [](double sum, double weight)
                                              { return sum + pow(weight, 2); }) /particles.size();
            adekf::viz::plotVector(Eigen::Matrix<double,1,1>(ess),"Ess",1000,"e");
            if (ess >= (double) cfg->lookup("ess") )
                return;
        }
        double u1 = distribution_(generator_);
        //probability step
        double one_by_n = 1. / particles.size();
        size_t particle_number = 0;
        double weight_sum = 0.;
        //temporary storage for the resampling
        adekf::aligned_vector<State> resampled_particles{particles.size()};
        std::vector<size_t> resampled_clusters(particles.size());
        std::vector<double> new_weights(weights.size());
        size_t index;
        
        //Iterate all particles
        for (index = 0; index < particles.size(); index++)
        {
            //Sum up weights up till now
            weight_sum += weights[index];

            //Add particle to resampled particles as long as the prob. sum in the resampled is smaller than prob. sum in original particles (stochastic resampling)
            while (u1 <= weight_sum)
            {
                u1 += one_by_n;//advance resampled prob. sum by 1/n
                new_weights[particle_number]=weights[index];
                resampled_particles[particle_number] = particles[index];
                resampled_clusters[particle_number++]= clusters[index];
            }

            if (u1 > 1.) //breaks if we would need a weight_sum  >1
            {
                
                break;
            }
        }
        //Check whether enough particles have been sampled. May happen due to  numeric errors. 
        if (particle_number <particles.size()){
            LOG(WARNING)<< "Did not sample enough particles: " << particles.size()-particle_number LOG_END;
            //Add missing particles from the back of the original
            while(particle_number < particles.size()){
                
                resampled_particles[particle_number] = particles.back();
                resampled_clusters[particle_number++]= clusters.back();
            }
        }
        //LOG_STREAM << "Breaking at: " << index LOG_END
        //write the resampled pdf to the input (inplace resampling)
        for (size_t i = 0; i < particles.size(); i++)
        {
            particles[i] = resampled_particles[i];
            clusters[i]=resampled_clusters[i];
            weights[i] = one_by_n;
        }
    }
    /**
     * @brief Shorthand to resample current pdf
     */
    void resample(bool use_ess)
    {
        resample(use_ess, weights_, particles_,clusters_);
    }

    
    /**
     * @brief Smoothes the old PDFs by Rauch Tung Striebel Smoothing.
     * 
     * This function takes considerable time and requires that the old PDFs where stored with storePDF.
     * 
     * It does not perform resampling of the old PDFs since this can deteoriate the old PDFs.
     * All old weights are normalized
     * 
     * @tparam DYNAMIC_MODEL  Class of the dynamic model
     * @tparam COV_TYPE  Type of the dynamic Covariance
     * @tparam Input Type of the dynamic model input
     * @tparam allocator Type of the input vector allocator
     * @param dynamicModel The dynamic Model of the system. Can only take 1 input. Use std::tuple to combine multiple.
     * @param covs The dynamic covariances for each time step. Have to be additive on the state. (predict can handle non additive as well)
     * @param inputs List of all dynamic Model inputs
     * @param input_start The index of the first input in the input list. I needed this, since I actually stored more inputs than I use. Can be ignored by most people.
     */
    template <class DYNAMIC_MODEL, class COV_TYPE, typename Input, typename allocator>
       void smooth(const DYNAMIC_MODEL &dynamicModel, const adekf::aligned_vector<COV_TYPE> &covs, const std::vector<Input, allocator> &inputs, size_t input_start=0){
        Eigen::Matrix<double,COV_TYPE::RowsAtCompileTime, 1> zero = zero.Zero();
        //iterate all old PDFs 
        for (int time_index = old_pdfs_.size() - 2; time_index >= 0; time_index--)
        {
            //Can calculate Gaussian probability between two particles
            mahalonobis_probability dynamic_probability{covs[time_index],false};
            // pdf at t+1
            auto n_plus_one_pdf = old_pdfs_[time_index + 1];
            // pdf at t
            auto n_pdf=old_pdfs_[time_index];
           
           //precompute dynamic probs for each particle of t to t+1
                Eigen::MatrixXd dyn_probs{weights_.size(),weights_.size()};
                for (size_t j = 0; j < weights_.size(); j++)
                {
                    for(size_t l=0; l < weights_.size();l++){
                        State prediction_l = n_pdf.particles_[l];
                        dynamicModel(prediction_l, zero, inputs[input_start+time_index]);
                        dyn_probs(l,j)=dynamic_probability(prediction_l, n_plus_one_pdf.particles_[j]);
                    }
                }
                
            std::vector<double> new_weights(n_pdf.weights_.size());
            //compute weight for each particle in t
             for (size_t i = 0; i < weights_.size(); i++)
            {
                

                //std::cout << prediction LOG_END;
                double new_weight = 0;
                for (size_t j = 0; j < weights_.size(); j++)
                {
                    // skip if the particle has 0 weight and thus no influence (breaks norm computation otherwise)
                    if (n_plus_one_pdf.weights_[j]==0)
                        continue; 
                    //Compute the normalization constant of the particle to avoid that information from before t influences t again.
                    double norm=0.;
                    for(size_t l=0; l < weights_.size();l++){
                        norm+=n_pdf.weights_[l]*dyn_probs(l,j);
                    }
                    //calculate the particle likelihood backpropagated from all time steps >t and all particles j
                    new_weight +=n_plus_one_pdf.weights_[j]* dyn_probs(i,j)/norm;
                }
                //Calculate smoothed weight as product forward pass likelihood and backward pass likelihood
                new_weights[i]=    n_pdf.weights_[i] * new_weight;
                DLOG(INFO) << time_index<<" "<<i<<" "<<new_weight LOG_END;
                
            }
            //Normalize all weights if weights did not deteriate due to some error
            //If error occurs the smoothing is interrupted and it is like its restarted at this time step
            if (normaliseWeights(new_weights)){
                for (size_t particle_index = 0; particle_index < weights_.size(); particle_index++){
                    old_pdfs_[time_index].weights_[particle_index]=new_weights[particle_index];
                } 
            }
        }
    }

    // Max iterations for getMean
    static inline size_t max_iterations = 100;
    // Break condition for getMean
    static inline double weighted_mean_epsilon = 1e-8;
    /**
     * @brief Computes the expected value of the PDF or of parts of it.
     * 
     * It can compute the expected value over parts of the PDF defined by the iterator begin and end
     * 
     * This function can handle vectors and manifolds. The approach is iterative.
     * For vectors computational overhead occurs.
     * 
     * Based on: https://www.mdpi.com/1424-8220/21/12/4164
     * 
     * @tparam Iter Class of the iterator
     * @param begin Iterator over indices to compute the mean of
     * @param end  End of begin
     * @param particles Particle states
     * @param weights Particle weights
     * @return State The expected value 
     */
    template<typename Iter>
    State getMean(Iter begin, Iter end, adekf::aligned_vector<State> & particles, std::vector<double> & weights)
    {

        State sum = particles[0];
        decltype(sum - sum) diff_sum = diff_sum.Zero();
        size_t iterations = 0;
        //iteratively computes the expected value
        do
        {
            //compute the expected difference between current mean and the pdf in the tangential space
            diff_sum = diff_sum.Zero();
            for (; begin < end; begin++)
            {
                diff_sum = diff_sum + weights[*begin] * (particles[*begin] - sum);
            }
            //adapt the mean until convergence or when to many iterations are used
            sum = sum + diff_sum; 
        } while (++iterations <= max_iterations && diff_sum.norm() > weighted_mean_epsilon);
        if (iterations > max_iterations)
            LOG(INFO) << "Warning: stopped due to excess of iterations" LOG_END;
         return sum;
    }

    /**
     * @brief Computes the sum of weights for parts of the weights
     * 
     * @tparam Iter Class of the iterator
     * @param begin Iterator that returns the desired indices
     * @param end End of begin
     * @param weights Particle weights
     * @return double Sum of weights
     */
    template<typename Iter>
    double getWeightSum(Iter begin, Iter end, std::vector<double> & weights)
    {

            double weight_sum = 0.;
            for (; begin < end; begin++)
            {
                weight_sum += weights[*begin];
            }
            return weight_sum;
    }

    /**
     * @brief Short hand to compute expected value of the whole current distribution
     * 
     * @return State Expected value of the current PDF
     */
    State getMean(){
        auto range = boost::irange(size_t(0), weights_.size());
        return getMean(range.begin(),range.end(),particles_,weights_);
    }
};

#endif //ZAVI_BOULDERN_CLUSTERPARTICLEFILTER_H
