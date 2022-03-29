/*
 * FileIMUPlugin.h
 *
 *  Created on: 16.01.2019
 *      Author: tomlucas
 */

#ifndef PLUGINS_FILEIMUPLUGIN_HPP_
#define PLUGINS_FILEIMUPLUGIN_HPP_
#define IMU_ROW_SIZE 14
#define OFFSET_TIL_START 5
#include "imu_plugin.hpp"

#include <iostream>
#include <fstream>

namespace zavi
::plugin {
	class FileIMUPlugin: public IMU_Plugin {
		typedef Eigen::Matrix<double,6,1> BIAS_MATRIX;
	public:
		bool reading_error=false;
		FileIMUPlugin(const std::string & filename, double skip_time=0, double end_time=-1.0, bool batchwise=false);
		virtual ~FileIMUPlugin();
		void run(double freq, bool (* running)());
		virtual bool can_batchwise() {return true;}

		/**
		 * @brief Loads all IMU data. 
		 * 
		 * Use when you want to use data multiple times 
		 * 
		 * @return int The amount of packages in the data
		 */
		size_t loadAll();

		/**
			 * @brief Get the time of the last data point
			 * 
			 * It is required to call @see loadAll before this.
			 * 
			 * @return double time of the last data point. -1 if loadAll is missing.
			 */
			double getLastTime();

		bool setTimeRange(double skip, double end);
		double getTimeAt(size_t index);


		bool readNextState();
		//Eigen::Vector3d acceleration_inc;
		size_t getCurrentRow() {
			return current_row;
		}
	protected:
		std::ifstream file;
		double skip;
		double end;
		std::vector<Eigen::Matrix<double, IMU_ROW_SIZE, 1>> rows;
		size_t current_row=0;
		const std::string name_of_file;
		const bool batchwise; // whether this one is in replay mode
		size_t added_states;
	};

}
/* namespace zavi::plugin */

#endif /* PLUGINS_FILEIMUPLUGIN_HPP_ */
