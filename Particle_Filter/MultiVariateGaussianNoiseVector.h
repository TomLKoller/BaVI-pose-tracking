//
// Created by tomlucas on 02.04.20.
//

#ifndef ADEKF_MULTIGAUSSIANNOISEVECTOR_H
#define ADEKF_MULTIGAUSSIANNOISEVECTOR_H

#include <Eigen/Core>
#include <random>
#include "ADEKFUtils.h"
#include <Eigen/Cholesky>
#include "probabilistics/mahalonobis_probability.hpp"
namespace adekf {

    /**
     * Class to create noise vectors with a gaussian random distribution.
     *
     * Use it like a usual eigen vector in computations. Call poll() to sample a new random vector.
     * @tparam Scalar type of the vector
     * @tparam size size of the vector
     */
    template<typename Scalar, int size>
    class MultiVariateGaussianNoiseVector : public Eigen::Matrix<Scalar, size, 1> {
        //Type of the vector
        using VectorType=Eigen::Matrix<Scalar, size, 1>;
        using Covariance=CovarianceOf<VectorType>;
        //generator for random numbers.
        static inline std::default_random_engine generator;
        //normal distribution to sample from
        std::normal_distribution<Scalar> distribution;
    public:
        //Mean and standard deviation of the gaussian distribution
        const VectorType mu;
        Covariance L;
        mahalonobis_probability<Covariance> probability_checker;

        /**
         * Initialises the vectors random distribution and assigns a random vector
         * @param mu mean of the gaussian
         * @param covariance standard deviation of the gaussian
         */
        MultiVariateGaussianNoiseVector(const VectorType & mu,const Covariance & covariance) : distribution(0, 1), mu(mu), L(covariance.llt().matrixL()), probability_checker(covariance) {
            this->poll();
        }
        /**
         * Samples a new gaussian random vector into this and also returns this
         */
        VectorType poll() {
            //Call move assign of the base class
            VectorType::operator=(mu+L*(this->unaryExpr([this](int) { return this->distribution(generator); })));
            return *this;
        }
        double probability(){
            return probability_checker(mu,*this);
        }


    };

    template<typename Derived>
    MultiVariateGaussianNoiseVector(const Derived & mu,const CovarianceOf<Derived> &cov ) ->MultiVariateGaussianNoiseVector<typename StateInfo<Derived>::ScalarType,DOFOf<Derived>>;
}


#endif //ADEKF_MULTIGAUSSIANNOISEVECTOR_H
