/*
 * state_plugin_base.hpp
 *
 *  Created on: 19.07.2018
 *      Author: tomlucas
 */

#ifndef PLUGINS_STATE_PLUGIN_BASE_HPP_
#define PLUGINS_STATE_PLUGIN_BASE_HPP_

#include <osg/Vec3d>
#include <osg/Quat>
namespace zavi{
namespace plugin{


/**
 * Converts an Eigen to an OSG Vector
 * @param eigen a eigen 3d vector
 * @return an osg 3d vector
 */
inline osg::Vec3d eigenToOSGVector(const Eigen::Vector3d & eigen) {
	return osg::Vec3d(eigen[0], eigen[1], eigen[2]);

}
/**
 * converts an ode quaternion to osg
 * @param ode_quat the ode quaternion
 * @return an osg quaternion
 */
inline osg::Quat eigenToOSGQuat(const Eigen::Quaterniond &eigen) {
	return osg::Quat(eigen.x(), eigen.y(), eigen.z(), eigen.w());     //<ODE uses  [w,x,y,z] and OSG [x,y,z,w]
}

/**
 * Base Class for all State plugins
 *
 *
 */
class StatePlugin{
public:
	/**
		 * Get the position as osg::Vec3D
		 * @return position of object
		 */
		virtual ::osg::Vec3d getPositionOSG()=0;

		/**
		 * Returns the rotation as osg::Quat
		 * @return the rotation
		 */
		virtual ::osg::Quat getRotationOSG()=0;

		virtual ~StatePlugin(){};

};

}//state_plugin
}//zavi




#endif /* PLUGINS_STATE_PLUGIN_BASE_HPP_ */
