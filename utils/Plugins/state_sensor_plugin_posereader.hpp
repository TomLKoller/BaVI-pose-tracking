/*
 * @file
 * state_plugin_estimator.hpp
 *
 *  Created on: 26.07.2018
 *      Author: tomlucas
 */

#ifndef PLUGINS_STATE_PLUGINSENSOR_POSEREADER_HPP_
#define PLUGINS_STATE_PLUGINSENSOR_POSEREADER_HPP_

#include <adekf_viz.h>

#include <Eigen/Core>
namespace zavi::plugin{

	template<class SensorType>
    class SensorPoseReader : public adekf::viz::PoseReader{
        SensorType * sensor;
    public:
        SensorPoseReader(SensorType * sensor):sensor(sensor){

        }
        virtual ~SensorPoseReader(){};

        Eigen::Vector3d getPosition() override {
            return sensor->getPosition();
        }

        Eigen::Quaterniond getOrientation() override {
            return sensor->getOrientation();
        }
    };


	template<class SensorType>
    class SensorPoseDelayedReader : public adekf::viz::PoseReader{
        SensorType * sensor;
        Eigen::Vector3d oldPos;
        Eigen::Vector3d showPos;
        Eigen::Quaterniond oldOrient;
        Eigen::Quaterniond showOrient;
    public:
        SensorPoseDelayedReader(SensorType * sensor):sensor(sensor),oldPos(sensor->getPosition()),showPos(oldPos),oldOrient(sensor->getOrientation()),showOrient(oldOrient){
            
        }
        virtual ~SensorPoseDelayedReader(){};

        Eigen::Vector3d getPosition() override {
            auto newPos=sensor->getPosition();
            if (oldPos != newPos){
                showPos=oldPos;
                oldPos=newPos;
            }
            return showPos;
        }

        Eigen::Quaterniond getOrientation() override {
            Eigen::Quaterniond newOrient= sensor->getOrientation();
            if(oldOrient.vec() !=newOrient.vec()){
                showOrient=oldOrient;
                oldOrient=newOrient;
            }
            return showOrient;
        }
    };

}//zavi::plugin




#endif /* PLUGINS_STATE_PLUGINSENSOR_POSEREADER_HPP_ */
