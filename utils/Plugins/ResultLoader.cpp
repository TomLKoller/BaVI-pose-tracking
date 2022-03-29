/*
 * ARTPlugin.cpp
 *
 *  Created on: 13.01.2020
 *      Author: tomlucas
 */

#include "ResultLoader.hpp"
#include "CSV_Reader.hpp"
#include "ZaviConfig.hpp"
#include <map>
#include <ADEKF/viz/adekf_viz.h>
namespace zavi ::plugin
{

	ResultLoader::ResultLoader(const char *filename, double skip_time, double end_time, bool batchwise) :  file(filename), skip(skip_time), end(end_time), name_of_file(filename), batchwise(batchwise)
	{
		if (end != -1)
			printf("Warning, End time not supported by ResultLoader");
		
		resetFile();

		LOG(INFO) <<"ResultLoader start time: " << nextTime << std::endl;
		//resetTime();
	}
	ResultLoader::~ResultLoader()
	{
	}
	
	bool ResultLoader::readNextState()
	{
		if (file.good())
		{
			std::string line;
			std::vector<std::string> splitted;
			getline(file, line);
			splitted = csv_reader::splitString(line, ',');
			if (splitted.size() > 1){
				return readLine(splitted);
			}
			else
			{
				LOG(WARNING) << "Wrong line format. Line:" << line;
				return false;
			}
		} 
		return false;
	}

	bool ResultLoader::readLine(std::vector<std::string> &splitted)
	{
		size_t loc_index = 0;
		size_t orient_index=4;
		if (splitted.size() > loc_index + 3)
		{
			position(0) = atof(splitted[loc_index + 1].c_str());
			position(1) = atof(splitted[loc_index + 2].c_str());
			position(2) = atof(splitted[loc_index + 3].c_str());
			nextTime = (atof(splitted[0].c_str()));
			imu_in_world=adekf::SO3d{atof(splitted[orient_index + 0].c_str()),atof(splitted[orient_index + 1].c_str()),atof(splitted[orient_index + 2].c_str()),atof(splitted[orient_index + 3].c_str())};
			//LOG(INFO) << "ART::Next time: " << nextTime LOG_END 
			return true;
		}
		else
		{
			LOG(WARNING) << "Failed to read a line at ART Plugin\n";
			return false;
		}
	}

	
	
	void ResultLoader::resetFile()
	{
		file.clear();
		file.seekg(0, file.beg);
		std::string skip_line{};
		getline(file, skip_line);//skip header line
		
		do {
			if(!readNextState()) {
				LOG(FATAL)<< "file stop before the given start time";
			}
		}while(nextTime<skip);
	}

	void ResultLoader::run(double freq, bool (*running)())
	{
		Timer loop(freq);
		double t_delta = 1. / (double)cfg->lookup("art_freq");
		newDataAvailable(t_delta);
		while (running())
		{
			if (readNextState())
			{
				Timer::setFileTime(oldTime);
				newDataAvailable(nextTime-oldTime);
				oldTime = nextTime;
				if (!batchwise)
					loop.wait();
			}
			else
			{
				batchReady();
				return;
			}
		}
	}

}
