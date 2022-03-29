#include "datasetFormat.hpp"
#include <fstream>
#include <iostream>
namespace bavi
{
    bool matchesRegex(const fs::directory_entry &entry, const std::regex &match)
{
    return (std::regex_match((std::string)entry.path().filename(), match));
}


bool isSearchedDirectory(const fs::directory_entry &entry,const std::regex &match)
{
    return entry.is_directory() && matchesRegex(entry, match);
}


std::list<fs::directory_entry> iterateFolder(const fs::path & folder,const std::regex &filter, bool (*test_func)(const fs::directory_entry &, const std::regex &)){
 std::list<fs::directory_entry> entries;
    for (const auto &entry : fs::directory_iterator(folder))
    {
        if ( test_func(entry,filter)){
            entries.push_back(entry);
        }
    }
    return entries;
}

std::list<fs::directory_entry> iterateParticipants(const fs::path & dataset_folder){
   return iterateFolder(dataset_folder,std::regex{"P\\d\\d"});
}

std::list<fs::directory_entry> iterateViews(const fs::path & participant_folder){
    return iterateFolder(participant_folder,std::regex("V\\d\\d(_\\d)?"));
}
std::list<fs::directory_entry> iterateTrials(const fs::path & view_folder){
    return iterateFolder(view_folder.string()
    +"/trials/",std::regex("T\\d+"));
}

json readJson(const fs::path &folder,const std::string & json_name){ 
        std::string file_name=folder.string() + "/"+json_name;
        std::ifstream json_file(file_name);   
        if (json_file.is_open()){
            json r_json;
            json_file >> r_json;
            return r_json;
        }
        else{
            std::cerr << "Json file " << file_name << " not found" <<std::endl;
            return NULL;
        }
        }

}//bavi 
