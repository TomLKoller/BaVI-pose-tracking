/*
 * ManifoldTypes.hpp
 *
 *  Created on: 19.12.2019
 *      Author: tomlucas
 */

#ifndef TYPES_MANIFOLDTYPES_HPP_
#define TYPES_MANIFOLDTYPES_HPP_

#include "ManifoldCreator.h"
#include "types/SO3.h"
namespace zavi::types{
	ADEKF_MANIFOLD(Pose3D,((adekf::SO3,orientation)),(3,position),(3,velocity))
	ADEKF_MANIFOLD(DoublePose3D,((Pose3D,first))((Pose3D, second)))
	ADEKF_MANIFOLD(XsensPose,((adekf::SO3,orientation)),(3,position),(3,velocity),(3,acc_bias),(3,ar_bias))
	ADEKF_MANIFOLD(XsensMarkerPose,((adekf::SO3,orientation)),(3,position),(3,velocity),(3,acc_bias),(3,ar_bias),(3,marker_offset))
	ADEKF_MANIFOLD(PoseVelo,,(3,position),(3,velocity))
}


#endif /* TYPES_MANIFOLDTYPES_HPP_ */
