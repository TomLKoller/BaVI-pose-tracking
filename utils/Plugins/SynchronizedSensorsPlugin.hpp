/*
 * SynchronizedSensorsPlugin.hpp
 *
 *  Created on: 17.01.2020
 *      Author: tomlucas
 */

#ifndef PLUGINS_SYNCHRONIZEDSENSORSPLUGIN_HPP_
#define PLUGINS_SYNCHRONIZEDSENSORSPLUGIN_HPP_

#include "sensor_plugin.hpp"
#include <vector>
namespace zavi
::plugin {

	class SynchronizedSensorPlugin :  public SensorPlugin {
	public:
		SynchronizedSensorPlugin(std::initializer_list<std::shared_ptr<SensorPlugin>> sensors, bool batchwise=true):sensors(sensors),batchwise(batchwise) {}

		/**
		 * Currently batch mode running  ( no sleep)
		 * @param freq
		 * @param running
		 */
		void run(double freq, bool (*running)()) {
			double current_time=0;
			auto loop=Timer(freq);
			while(running() ) {
				for(auto sensor : sensors) {
					// update current sensor
					if(sensor->nextTime == current_time) {
						Timer::setFileTime(current_time);
						sensor->newDataAvailable(sensor->nextTime-sensor->oldTime);
						newDataAvailable(sensor->nextTime-sensor->oldTime);
						sensor->oldTime=current_time;
						sensor->readNextState();
					}

				}
				//find next sensor update time
				double min_time=std::numeric_limits<double>::max();

				for(auto sensor: sensors){
					if(sensor->nextTime > current_time and sensor->nextTime < min_time)
						min_time=sensor->nextTime;
				}
				if(min_time > current_time)
					current_time=min_time;
				else // no more valid data
					break;
				if(!batchwise)
				loop.wait();

			}
			batchReady();
			LOG(INFO) << "Finished Synchronized Sensors\n";
			zavi::printf("Finished Sensors");

		}

	private:
		std::vector<std::shared_ptr<SensorPlugin>> sensors;
		bool batchwise;
	};

}

#endif /* PLUGINS_SYNCHRONIZEDSENSORSPLUGIN_HPP_ */
