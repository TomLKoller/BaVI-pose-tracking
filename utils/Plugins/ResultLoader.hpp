/*
 * ARTPlugin.hpp
 *
 *  Created on: 13.01.2020
 *      Author: tomlucas
 */

#ifndef PLUGINS_GTLOADER_HPP_
#define PLUGINS_GTLOADER_HPP_


#include "sensor_plugin.hpp"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include  <ADEKF/types/SO3.h>
#include <set>
#include <fstream>
namespace zavi::plugin{
	class ResultLoader : public SensorPlugin{
	public:
		ResultLoader(const char * filename,double skip_time=0, double end_time=-1.0, bool batchwise=false);
			virtual ~ResultLoader();
			void run(double freq, bool (* running)());
			virtual bool can_batchwise() {return true;}
			bool readNextState();

			
			void exportData(const char * export_file);

			Eigen::Vector3d getPosition(){
				return position;
			}
			adekf::SO3<double> getOrientation(){
				return imu_in_world;

			}


			void resetTime(){
				skip=nextTime;
				nextTime=0;
				oldTime=0;
			}

			void resetFile();

		protected:

			Eigen::Vector3d position;
			adekf::SO3<double> imu_in_world;
			std::ifstream file;
			double skip;
			double end;
			const char * name_of_file;
			const bool batchwise; // whether this one is in replay mode

			/**
			 * Reads a line to the attributes of the plugin
			 * @return success of reading the data
			 */
			bool readLine(std::vector<std::string> &splitted);


	};

}


#endif /* PLUGINS_GTLOADER_HPP_ */
