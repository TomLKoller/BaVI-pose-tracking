/**
 * @file gripDetector.h 
 * @author Tom Koller tkoller@uni-bremen.de
 * @brief This file contains the grip Detector for Bouldering
 * @version 0.1
 * @date 2022-03-17
 * 
 * 
 */
#include <vector>
#include "Plugins/FileIMUPlugin.hpp"
#include "ZaVI_Utils.hpp"


/**
 * @brief Detects grips in the IMU data.
 * 
 * Calculates angular rate std over moving windows and performs a thresholding.
 * The windows are centered around the current time index, i.e. [k-width/2,k+width/2]
 * 
 * The first and last width/2 time steps are always false since no window can be found
 * 
 *Based on: C. Ladha, N. Y. Hammerla, P. Olivier, and T. Plötz, “ClimbAX: skill
 *assessment for climbing enthusiasts,” in Proceedings of the 2013 ACM
 *international joint conference on Pervasive and ubiquitous computing,
 *New York, NY, USA, Sep. 2013, pp. 235–244. https://doi.org/10.1145/2493432.2493492
 * 
 * @param imu IMU plugin with already loaded data and correct time range
 * @param thresh treshhold to detect low acceleration phases
 * @param width Width of the moving window to accumulate acceleration 
 * @param center_jump Offset between windows
 * @return std::vector<bool> A vector with true for grip and false for no grip
 */
std::vector<bool> gripDetector(std::shared_ptr<zavi::plugin::FileIMUPlugin> imu, double thresh, size_t width=28, size_t center_jump=1)
{
    //Read World accelerations and angular rates based on xsens orientation
     adekf::aligned_vector<Eigen::Vector3d> accelerations;
     adekf::aligned_vector<Eigen::Vector3d> angular_rates;
    while (imu->readNextState())
    {
        accelerations.push_back(imu->getOrientation() * imu->getAcceleration() - zavi::plugin::IMU_Plugin::gravity);
        angular_rates.push_back(imu->getOrientation() * imu->getAngularVelocity());
    }

    std::vector<bool> is_grip;
    for (size_t j = 0; j < width / 2; j++) // start with false since no window possible
        is_grip.push_back(0);

    //Find all grips
    for (size_t i = width / 2; i < accelerations.size() - width / 2; i += center_jump)
    {
        //Calculate mean AR of window
        Eigen::Vector3d meanA = meanA.Zero(), meanAR = meanA.Zero();
        for (size_t j = i - width / 2; j < i + width / 2; j++)
        {
            //meanA += accelerations[j];
            meanAR += angular_rates[j];
        }
        //meanA /= width;
        meanAR /= width;
        //Calculate standard deviations
        Eigen::Array3d stdA = stdA.Zero(), stdAR = stdAR.Zero();
        for (size_t j = i - width / 2; j < i + width / 2; j++)
        {
            //stdA += (accelerations[j] - meanA).array().square();
            stdAR += (angular_rates[j] - meanAR).array().square();
        }
        //stdA /= width;
        stdAR /= width;
        //stdA = stdA.sqrt();
        stdAR = stdAR.sqrt();
        is_grip.push_back(stdAR.sum() < thresh);

    }

    for (size_t j = 0; j < width / 2; j++) // end since no window possible
        is_grip.push_back(0);
    return is_grip;
}