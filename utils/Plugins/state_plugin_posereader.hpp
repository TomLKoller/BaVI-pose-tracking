/*
 * @file
 * state_plugin_estimator.hpp
 *
 *  Created on: 26.07.2018
 *      Author: tomlucas
 */

#ifndef PLUGINS_STATE_PLUGIN_POSEREADER_HPP_
#define PLUGINS_STATE_PLUGIN_POSEREADER_HPP_

#include "state_plugin_base.hpp"
#include "PoseReader.hpp"
#include <zavi/utils/OSG_Utils.hpp>

#include <Eigen/Core>
namespace zavi::plugin{


	class StatePluginEstimator : public StatePlugin{
public:

				StatePluginEstimator(std::shared_ptr<plugin::PoseReader> estimator):estimator(estimator){

				}
				inline ::osg::Vec3d getPositionOSG(){
						return zavi::osg_viz::eigenToOSGVector(  estimator->getEstimatedPosition());
				}
				inline Eigen::Matrix3d getPositionError(){
						return estimator->getEstimatedPositionError();
				}


				inline ::osg::Quat getRotationOSG(){
					return zavi::osg_viz::eigenToOSGQuat(estimator->getEstimatedOrientation());
				}

private:
				std::shared_ptr<plugin::PoseReader> estimator;
	};

}//zavi::plugin




#endif /* PLUGINS_STATE_PLUGIN_POSEREADER_HPP_ */
