/*
 * ARTPlugin.cpp
 *
 *  Created on: 13.01.2020
 *      Author: tomlucas
 */

#include "ARTPlugin.hpp"
#include "CSV_Reader.hpp"
#include "ZaviConfig.hpp"
#include <map>
#include <ADEKF/viz/adekf_viz.h>
namespace zavi ::plugin
{

	ARTPlugin::ARTPlugin(const char *filename, std::set<int> markerIDs, Eigen::Quaterniond art_in_world_rotation, const Eigen::Vector3d &art_in_world_translation, double skip_time, double end_time, bool batchwise) : markerIDs(markerIDs), file(filename), art_in_world_rotation(art_in_world_rotation), art_in_world_translation(art_in_world_translation), marker_id(-1), skip(0), end(end_time), name_of_file(filename), batchwise(batchwise)
	{
		if (skip_time != 0)
			printf("WARNING, Skip time not supported by ART");
		if (end != -1)
			printf("Warning, End time not supported by ART");
		resetFile();

		LOG(INFO) << "Art start time: " << nextTime << std::endl;
		//resetTime();
		orientation = orientation.Identity();
	}
	ARTPlugin::~ARTPlugin()
	{
	}
	int readMarkerID(const std::vector<std::string> &splitted)
	{
		return atoi(splitted[1].c_str());
	}

	bool ARTPlugin::readNextState()
	{

		std::string line;
		std::vector<std::string> splitted;
		Eigen::Matrix<double, ART_ROW_SIZE, 1> row;
		do
		{
			if (current_row < rows.size())
			{
				row = rows[current_row];
				current_row++;
			}
			//Else load new row from file

			else if (!zavi::csv_reader::readLineToEigen(file, row))
			{
				return false;
			}

			marker_id = row(1);
		} while ((!markerIDs.empty() && markerIDs.find(marker_id) == markerIDs.end()) && file.good()); //Skip all lines with unwanted markers
		//Check if at end of file without finding a wanted marker
		if (!markerIDs.empty() && markerIDs.find(marker_id) == markerIDs.end())
			return false;

		return readLine(row);

		return false;
	}

	size_t ARTPlugin::loadAll()
	{
		LOG(INFO) << "Loading data of Art" << std::endl;
		if (skip != 0. || end != -1)
		{
			LOG(WARNING) << "loadAll is not intended for use with skip(" << skip << ") and end(" << end << ")" << std::endl;
		}
		Eigen::Matrix<double, ART_ROW_SIZE, 1> row;
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

	double ARTPlugin::getLastTime()
	{
		if (rows.size() == 0)
		{
			LOG(ERROR) << "Tried to retreive last time from ART without loading all data first" << std::endl;
			return -1;
		}
		return rows.back()(0) - skip;
	}

	bool ARTPlugin::readLine(const Eigen::Matrix<double, ART_ROW_SIZE, 1> &row)
	{
		size_t loc_index = 1;

		position = row.segment<3>(loc_index + 1);
		nextTime = row(0) - skip;
		if (!(end == -1. or nextTime <= end))
		{
			return false;
		}
		position = art_in_world_rotation * position - art_in_world_translation;
		//LOG(INFO) << "ART::Next time: " << nextTime LOG_END
		return true;
	}

	void ARTPlugin::showBubbles(double size)
	{
		std::map<int, const char *> color_map;
		std::map<size_t, std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>> positions;
		resetFile();
		do
		{
			if (color_map.find(marker_id) == color_map.end())
				color_map.emplace(marker_id, adekf::viz::colorpalette[color_map.size() % 7]);

			positions[marker_id].push_back(getPosition());
		} while (readNextState());
		for (auto tuple : positions)
		{
			//adekf::viz::PoseRenderer::displayPath(tuple.second, color_map[tuple.first]);
			adekf::viz::PoseRenderer::displayPoints(tuple.second, color_map[tuple.first], size);
		}
		resetFile();
		printf("Finished adding bubbles");
	}

	void ARTPlugin::exportData(const char *export_file)
	{
		std::ofstream ex_file(export_file);
		if (!ex_file.is_open()){
			std::cout << "Could not open File" <<std::endl;
			return;
		}
		resetFile();
		ex_file << "time in s, marker_id, x in m, y in m, z in m" << std::endl;
		do
		{
			ex_file << nextTime << "," << marker_id << "," << position.x() << "," << position.y() << "," << position.z() << std::endl;
		} while (readNextState());
		ex_file.close();
	}

	void ARTPlugin::resetFile()
	{
		file.clear();
		file.seekg(0, file.beg);
		std::string skip_line{};
		getline(file, skip_line);
		readNextState();
		skip = 0;
	}

	void ARTPlugin::setTimeRange(double skip, double end)
	{
		size_t i = 0;
		while (nextTime < skip && readNextState())
		{
			i++;
		}
		if (!file.good())
			LOG(ERROR) << "ARTPlugin: Requested time range unavailable" << std::endl;
		//this->skip=nextTime;
		//Enable end =-1 to use the rest of the file
		this->end = end;
		//nextTime=skip;
		//set start values
		readNextState();
		oldTime = skip;
	}

	void ARTPlugin::run(double freq, bool (*running)())
	{
		Timer loop(freq);
		double t_delta = 1. / (double)cfg->lookup("art_freq");
		newDataAvailable(t_delta);
		while (running())
		{
			if (readNextState())
			{
				Timer::setFileTime(oldTime);
				oldTime = nextTime;
				newDataAvailable(t_delta);
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
