//
// Setup to store grip detection and to check trials for outliers
// Compile this instead of realDataSmoothing.cpp to check for outliers
// Created by Tom Koller on 03.11.20.
//
#define SHOW_VIZ
#include <filesystem>
#include <nlohmann/json.hpp>
#include "ZaviConfig.hpp"
#include "datasetFormat.hpp"
#include "ZaVI_Utils.hpp"
#include "Plugins/FileIMUPlugin.hpp"
#include "gripDetector.h"
namespace fs = std::filesystem;
using json = nlohmann::json;
using namespace Eigen;

/**
 * @brief Checks whether a trial has outliers and exports the detected grips
 *
 * @param imu The IMU plugin of the trial
 * @param trial_data The trial_data.json data
 * @param grip_output_path The output path for the grip detector
 * @return true If has IMU outlier
 * @return false Else
 */
bool hasOutlier(std::shared_ptr<zavi::plugin::FileIMUPlugin> imu, const bavi::json &trial_data, std::string &grip_output_path)
{
    // read times from trial json
    double start_time = trial_data["start_time"];
    double end_time = trial_data["end_time"];
    LOG(INFO) << "Trial starts at: " << start_time LOG_END;
    imu->setTimeRange(start_time, end_time);
    size_t imu_current_row = imu->getCurrentRow();
    auto is_grip = gripDetector(imu, cfg->lookup("std_thresh"));
    LOG(INFO) << "Grip vector size: " << is_grip.size() LOG_END;

    // Find last grip
    size_t last_grip = 0;
    for (size_t index = 0; index < is_grip.size(); index++)
    {
        if (is_grip[index])
        {
            last_grip = index;
        }
    }
    // Get trial end at last grip
    end_time = imu->getTimeAt(imu_current_row + last_grip);

    // Write all grips to output
    std::ofstream grip_output{grip_output_path};
    grip_output << "Time in s, grip" << std::endl;

    for (size_t index = 0; index < last_grip; index++)
    {
        grip_output << imu->getTimeAt(imu_current_row + index) << "," << (is_grip[index] ? 1 : 0) << std::endl;
    }
    grip_output.flush();
    grip_output.close();

    // Check for all readings if outlier occured
    imu->setTimeRange(start_time, end_time); // imu data is used in gripDetector
    while (imu->readNextState())
    {
        if (imu->reading_error)
            return true;
    }
    return false;
}
/**
 * @brief Traverses the dataset and finds outlier trials
 *
 * @param argc
 * @param argv
 * @return int
 */
int main(int argc, char **argv)
{
    // Check argument format
    if (argc != 7)
    {
        std::cout << "Passed: ";
        for (size_t i = 0; i < argc; i++)
            std::cout << argv[i] << " ";
        std::cout << std::endl;
        std::cout << "Expected Format: Zavi-Particle_Filter_Smooth <dataset_directory>  <result_directory> <smooth_result_directory> <transition_result_directory> <map_directory> <Config_folder>" << std::endl;
        return -1;
    }
    printf("Starting up\n");
    initConfig(argv[6], "main_config.cfg");
    google::InitGoogleLogging(argv[0]);
    google::SetStderrLogging(cfg->lookup("log_level"));
    zavi::Timer::initTime(cfg->lookup("time_lapse"));
    // Check if only specific trials are wanted
    std::string target_trial = cfg->lookup("target_trial");
    std::string target_sensor = cfg->lookup("target_sensor");

    std::filesystem::path dataset_path{argv[1]};
    //Traverse all participants
    for (const auto &participant_entry : bavi::iterateParticipants(dataset_path))
    {
        std::filesystem::path participant_folder{participant_entry.path()};
        participant_folder = participant_folder.lexically_normal();
        std::cout << participant_folder << std::endl;

        //Copy participant_data.json to result folder
        std::filesystem::path participant_data{participant_folder.string() + "/participant_data.json"};
        participant_data = participant_data.lexically_normal();
        for (char *folder : {argv[2], argv[3], argv[4]})
        {
            std::filesystem::path participant_data_output = std::string(folder) + (participant_data).lexically_relative(dataset_path).string();
            fs::create_directories(participant_data_output.parent_path());
            std::filesystem::copy_file(participant_data, participant_data_output, std::filesystem::copy_options::update_existing);
        }

       
        // Iterate all viewpoints
        for (const auto &view_folder : bavi::iterateViews(participant_folder))
        {
            // Iterate all trials
            for (const auto &trial_folder : bavi::iterateTrials(view_folder.path()))
            {
                //Check if trial is wanted
                if (target_trial != "_" && trial_folder.path().string().find(target_trial) == std::string::npos)
                    continue; // Skip if not the target trial
                std::cout << trial_folder LOG_END;

                // Copy trial_data.json
                std::filesystem::path trial_data_file = trial_folder.path().string() + "/trial_data.json";
                trial_data_file = trial_data_file.lexically_normal();
                for (char *folder : {argv[2], argv[3], argv[4]})
                {

                    std::filesystem::path trial_data_output = std::string(folder) + std::filesystem::path(trial_data_file).lexically_relative(dataset_path).string();
                    fs::create_directories(trial_data_output.parent_path());
                    std::filesystem::copy_file(trial_data_file, trial_data_output, std::filesystem::copy_options::overwrite_existing);
                }
                //read trial data
                auto trial_data = bavi::readJson(trial_folder, "trial_data.json");
                //Check if synched video is available, otherwise this trial does not have start time end time etc. 
                if (trial_data["video_synched"] == 0)
                {
                    LOG(INFO) << "Skipped due to missing video (no start_time available)" LOG_END;
                    continue;
                }

                // Find all imu_files
                for (const auto &imu_file : bavi::iterateFolder(trial_folder, std::regex{"imu_\\w+.csv"}, bavi::matchesRegex))
                {
                    // Skip if not the target sensor
                    if (target_sensor != "_" && imu_file.path().string().find(target_sensor) == std::string::npos)
                        continue; 
                    auto imu = std::make_shared<zavi::plugin::FileIMUPlugin>(imu_file.path().string(), 0, -1, false);
                    imu->loadAll();
                    //create output name 
                    std::regex ex{"(.+?)(imu_)(\\w+)(.csv)"};
                    std::smatch match;
                    std::string temp = imu_file.path().string();
                    std::string output_name = std::regex_replace(imu_file.path().string(), ex, "$1grip/grip_$3.csv");
                    // std::cout << output_name LOG_END;
                    fs::create_directories(std::filesystem::path{output_name}.parent_path());
                    //Extract the sensor name
                    if (std::regex_search(temp, match, ex))
                    {
                        std::string sensor_name = match[3];
                        // std::cout << sensor_name LOG_END;
                        int result = hasOutlier(imu, trial_data, output_name) ? 1 : 0;
                        // std::cout << result LOG_END;
                        //Set outlier flag in trial data
                        trial_data["has_outliers"][sensor_name] = result;
                    }
                    // return;
                }
                //write new trial data
                std::string trial_data_input = trial_folder.path().string() + "/trial_data.json";
                std::ofstream out_stream(trial_data_input);
                out_stream << std::setw(4) << trial_data;
                out_stream.close();
            }
        }
    }
    return 0;
}