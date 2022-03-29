/*
 * FileIMUPlugin.cpp
 *
 *  Created on: 16.01.2019
 *      Author: tomlucas
 */

#include "FileIMUPlugin.hpp"
#include "CSV_Reader.hpp"
#include "ZaVI_Utils.hpp"
#include "ZaviConfig.hpp"
#include <glog/logging.h>
namespace zavi ::plugin
{

	FileIMUPlugin::FileIMUPlugin(const std::string &filename, double skip_time, double end_time, bool batchwise) : file(filename), skip(skip_time), name_of_file(filename), batchwise(batchwise), added_states(0)
	{
		Timer::activateFileMode();
		acceleration_cov = Eigen::Matrix3d::Identity() * 0.0004;
		orientation_cov = Eigen::Matrix3d::Zero();
		orientation_cov(0, 0) = pow(0.75 * M_PI / 180, 2); //from website
		orientation_cov(1, 1) = pow(0.75 * M_PI / 180, 2); //from website
		orientation_cov(2, 2) = pow(1.5 * M_PI / 180, 2);  //from website
		angular_velocity_cov = Eigen::Matrix3d::Identity() * pow(0.1 * M_PI / 180., 2);
		acceleration_bias_cov = acceleration_bias_cov.Identity() * pow(0.1 * 0.91e-3, 2); // bias stability from whitepaper
		acceleration_bias_corr = BIAS_MATRIX(std::initializer_list<double>({1000, 2000, 1000, 2000, 1000, 2000}).begin());
		angular_rate_bias_cov = angular_rate_bias_cov.Identity() * pow(10. / 180. * M_PI / 3600., 2);
		angular_rate_bias_corr = Eigen::Vector3d(2000, 2000, 2000);
		// TODO Auto-generated constructor stub
		std::string line;

		std::getline(file, line); // skip first line
		Eigen::Matrix<double, IMU_ROW_SIZE, 1> row = row.Zero();
		do
		{
			if (!zavi::csv_reader::readLineToEigen(file, row))
			{
				LOG(FATAL) << "file stop before the given start time";
			}
			current_row++;
		} while (row(0) < skip_time);
		current_row--;
		this->skip = row(0);
		LOG(INFO) << "Skip is " << skip;
		angular_velocity = row.block(4, 0, 3, 1);
		imu_in_world = Eigen::Quaterniond(row(10), row(11), row(12), row(13)); //Caution, Eigen internally stores as xyzw wherefore mapping does not work.
		acceleration = row.block(1, 0, 3, 1);
		magnetometer = row.block(7, 0, 3, 1);
		nextTime = this->skip;
		oldTime = this->skip - 1. / (double)cfg->lookup("imu_freq");
		end = end_time;
	}

	FileIMUPlugin::~FileIMUPlugin()
	{
		DLOG(INFO) << "File imu got destroyed";
	}

	bool FileIMUPlugin::readNextState()
	{
		Eigen::Matrix<double, IMU_ROW_SIZE, 1> row;
		//LOG(INFO) << "Current row: " << current_row <<std::endl;
		//if already read load from rows
		if (current_row < rows.size())
		{
			row = rows[current_row];
		}
		//Else load new row from file
		else if ((end == -1. or nextTime <= end) and file.good())
		{
			if (!zavi::csv_reader::readLineToEigen(file, row))
			{
				return false;
			}
		}
		else
		{
			return false; // if both not viable
		}
		reading_error=false;
		//LOG(INFO) << "IMU row:" << row.transpose() <<std::endl;
		nextTime = row(0);
		if (nextTime < 0. || (end != -1 && nextTime > end))
		{
			return false; //out of skip-end region
		}
		acceleration = row.block(1, 0, 3, 1);
		angular_velocity = row.block(4, 0, 3, 1);
		if (acceleration.norm() > 160.)
		{
			LOG(INFO) << row;
			LOG(WARNING) << "Found acceleration: " << acceleration.norm() << " It is over sensor limit. Found at: " << row(0);
			reading_error=true;
		}
		if (angular_velocity.norm() > 2000. / 180. * M_PI)
		{
			LOG(INFO) << row;
			LOG(WARNING) << "Found angular velocity: " << angular_velocity.norm() << " It is over sensor limit. Found at: " << row(0);
			reading_error=true;
		}
		magnetometer = row.block(7, 0, 3, 1);
		if (magnetometer.norm() > 1.9)
		{
			LOG(INFO) << row;
			LOG(WARNING) << "Found magnetic field of: " << magnetometer.norm() << " It is over sensor limit. Found at: " << row(0);
		}
		imu_in_world = Eigen::Quaterniond(row(10), row(11), row(12), row(13)); //Caution, Eigen internally stores as xyzw wherefore mapping does not work.
		current_row++;
		return true;
	}
	size_t FileIMUPlugin::loadAll()
	{
		LOG(INFO) << "Loading data of IMU" << std::endl;
		if (skip != 0. || end != -1)
		{
			LOG(WARNING) << "loadAll is not intended for use with skip(" << skip << ") and end(" << end << ")" << std::endl;
		}
		Eigen::Matrix<double, IMU_ROW_SIZE, 1> row;
		file.clear();
		file.seekg(0, file.beg);
		rows.clear();
		//skip first line
		zavi::csv_reader::readLineToEigen(file, row);
		while (file.good())
		{
			if (!zavi::csv_reader::readLineToEigen(file, row))
			{

				LOG(INFO) << "Loaded " << rows.size() << " lines." << std::endl;
				break;
			}
			rows.push_back(row);
		}
		current_row = 0;
		readNextState();

		return rows.size();
	}

		double FileIMUPlugin::getLastTime()
	{
		if (rows.size() == 0)
		{
			LOG(ERROR) << "Tried to retreive last time from ART without loading all data first" << std::endl;
			return -1;
		}
		if (end ==-1)
			return rows.back()(0);
		else 
			return end;
	}
	double FileIMUPlugin::getTimeAt(size_t index){
		return rows[index](0);
	}


	bool FileIMUPlugin::setTimeRange(double skip, double end)
	{
		size_t i = 0;
		while (rows[i](0) < skip)
		{
			i++;
		}
		if (i >= rows.size()){
			LOG(ERROR) << "FileIMU: Requested time range unavailable" << std::endl;
			return false;
		}
		this->skip = rows[i](0);
		this->end = end;
		current_row = i;
		nextTime = skip;
		//set start values
		readNextState();
		if (i > 0)
		{
			oldTime = rows[i - 1](0);
		}
		else
		{
			oldTime = skip - 1. / (double)cfg->lookup("imu_freq");
		}
		return true;
	}
	void FileIMUPlugin::run(double freq, bool (*running)())
	{
		Timer loop(freq);
		//create noises
		oldTime = skip;
		double t_delta = 1. / (double)cfg->lookup("imu_freq");
		// fixes that i read the first value in the constructor
		newDataAvailable(t_delta);
		size_t calls = 0;
		while (running())
		{
			if (readNextState())
			{

				calls++;
				double t_diff = nextTime - oldTime - t_delta;
				while (round(t_diff / t_delta) > 0)
				{
					data_valid = false;
					Timer::setFileTime(nextTime - t_diff);
					newDataAvailable(nextTime - oldTime - t_diff);
					t_diff -= t_delta;
					added_states++;
				}
				data_valid = true;
				Timer::setFileTime(nextTime);
				newDataAvailable(nextTime - oldTime);
				oldTime = nextTime;
				if (!batchwise)
					loop.wait();
			}
			else
			{
				batchReady();
				LOG(INFO) << "Finished file " << this->name_of_file << " at " << oldTime << "s";
				LOG(INFO) << "Read  " << calls << " valid states from file";
				LOG(INFO) << "Added " << added_states << " states that have no valid input" << std::endl;
				return;
			}
		}
	}
}
//namespace zavi::plugin
