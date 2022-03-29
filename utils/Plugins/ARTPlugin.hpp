/*
 * ARTPlugin.hpp
 *
 *  Created on: 13.01.2020
 *      Author: tomlucas
 */

#ifndef PLUGINS_ARTPLUGIN_HPP_
#define PLUGINS_ARTPLUGIN_HPP_
#define ART_ROW_SIZE 5

#include "sensor_plugin.hpp"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include  <ADEKF/types/SO3.h>
#include <set>
#include <fstream>
namespace zavi::plugin{
	class ARTPlugin : public SensorPlugin{
	public:
		ARTPlugin(const char * filename,std::set<int> markerIDs,Eigen::Quaterniond art_in_world=Eigen::Quaterniond::Identity(),const Eigen::Vector3d &art_in_world_translation=Eigen::Vector3d::Zero(), double skip_time=0, double end_time=-1.0, bool batchwise=false);
			virtual ~ARTPlugin();
			void run(double freq, bool (* running)());
			virtual bool can_batchwise() {return true;}
			bool readNextState();

			void showBubbles(double size);

			void exportData(const char * export_file);

			Eigen::Vector3d getPosition(){
				return position;
			}
			int getMarker(){ return marker_id;}
			adekf::SO3<double> getOrientation(){
				return orientation;

			}


		/**
		 * @brief Loads all Art data. 
		 * 
		 * Use when you want to use data multiple times 
		 * 
		 * @return int The amount of packages in the data
		 */
		size_t loadAll();

			void resetTime(){
				skip=nextTime;
				nextTime=0;
				oldTime=0;
			}

			void setTimeRange(double start, double end=1);

			void resetFile();
			/**
			 * @brief Get the time of the last data point
			 * 
			 * It is required to call @see loadAll before this.
			 * 
			 * @return double time of the last data point. -1 if loadAll is missing.
			 */
			double getLastTime();

			//!All IDs belonging to this object are saved here
			std::set<int> markerIDs;
		protected:

			Eigen::Vector3d position;
			adekf::SO3<double> orientation;
			std::ifstream file;
			Eigen::Quaterniond art_in_world_rotation;
			Eigen::Vector3d art_in_world_translation;
			int marker_id;
			double skip;
			double end;
			const char * name_of_file;
			const bool batchwise; // whether this one is in replay mode
			std::vector<Eigen::Matrix<double, ART_ROW_SIZE, 1>> rows;
			size_t current_row=0;
		

			/**
			 * Reads a line to the attributes of the plugin
			 * @return success of reading the data
			 */
			bool readLine(const Eigen::Matrix<double, ART_ROW_SIZE, 1> & row);


	};

}


#endif /* PLUGINS_ARTPLUGIN_HPP_ */
