#pragma once
/**
 * Struct to calculate the probability on a multivariate gaussian
 * All results are logarithmically
 */
template<typename STIFF>
struct mahalonobis_probability {
    //The inversed covariance
    typename STIFF::PlainObject cov_inv;
    //Whether the output is in log probability or not (defaults to true)
    bool log_output;
    //The factor in front of the exponential (logarithmic)
    double factor;
    
    /**
     * Calculates the inverse covariance and the factor of the gaussian
     * @param cov  The covariance of the Gaussian distribution
     */
    mahalonobis_probability(const Eigen::MatrixBase<STIFF> &cov, bool log_output=true) : cov_inv(cov.inverse()), log_output(log_output) {
        factor = log(1. / (pow(2 * M_PI, STIFF::RowsAtCompileTime / 2.) * sqrt(abs(cov.determinant()) )));
    }

    /**
     * Calculates the logarithmic probability of  (a-b)
     * @tparam Derived The type of a and b
     * @param a minuend
     * @param b subtrahend
     * @return log(N((a-b),cov))
     */
    template<typename DerivedA, typename DerivedB>
    auto operator()(const DerivedA &a, const DerivedB &b) {
        auto diff = (a - b).eval();
        auto result= factor - 0.5 * (diff.transpose() * cov_inv * diff)(0);
        if (!log_output)
            return exp(result);
        else
            return result;
    }
};


/**
 * Struct to calculate the probability on a multivariate gaussian
 * All results are logarithmically
 */
template<typename STIFF>
struct mahalonobis_probability_with_set {
    //The inversed covariance
    typename STIFF::PlainObject cov_inv;
    double circle_region;
    //Whether the output is in log probability or not
    bool log_output;
     //The factor in front of the exponential (logarithmic)
    double factor;
    double min_return;
    
    /**
     * Calculates the inverse covariance and the factor of the gaussian
     * @param cov  The covariance of the Gaussian distribution
     */
    mahalonobis_probability_with_set(const Eigen::MatrixBase<STIFF> &cov, const double circle_region, const double num_stds, bool log_output=true) : cov_inv(cov.inverse()), circle_region(circle_region), log_output(log_output) {
        factor = log(1. / (pow(2 * M_PI, STIFF::RowsAtCompileTime / 2.) * sqrt(abs(cov.determinant()) )));
        min_return=factor-0.5*pow(num_stds,2); //  we do not distinguish any further than the num_stds
    }

    /**
     * Calculates the logarithmic probability of  (a-b)
     * @tparam Derived The type of a and b
     * @param a minuend
     * @param b subtrahend
     * @return log(N((a-b),cov))
     */
    template<typename DerivedA, typename DerivedB>
    double operator()(const DerivedA &a, const DerivedB &b) {
        auto diff = (a - b).eval();
        auto diff_norm=diff.norm();
        if (diff_norm< circle_region)
             return factor;
        diff-=diff/diff_norm*circle_region;
        auto result= std::max(factor - 0.5 * (diff.transpose() * cov_inv * diff)(0),min_return);
        if (!log_output)
            return exp(result);
        else
            return result;
    }
};