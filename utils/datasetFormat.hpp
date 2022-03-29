/*
 * Callbacker.hpp
 *
 *  Created on: 25.06.2021
 *      Author: tomlucas
 */

#ifndef DATASETFORMAT_HPP_
#define DATASETFORMAT_HPP_

#include <filesystem>
#include <nlohmann/json.hpp>
#include <regex>
#include <list>

namespace bavi{
namespace fs = std::filesystem;

using json = ::nlohmann::json;

bool matchesRegex(const fs::directory_entry &entry,const std::regex &match);

bool isSearchedDirectory(const fs::directory_entry &entry,const std::regex &match);
std::list<fs::directory_entry> iterateFolder(const fs::path & folder,const std::regex &filter, bool (*test_func)(const fs::directory_entry &, const std::regex &)=isSearchedDirectory);
std::list<fs::directory_entry> iterateParticipants(const fs::path & dataset_folder);

std::list<fs::directory_entry> iterateViews(const fs::path & participant_folder);
std::list<fs::directory_entry> iterateTrials(const fs::path & view_folder);

json readJson(const fs::path &folder,const std::string & json_name);
}//bavi

#endif /* DATASETFORMAT_HPP_ */
