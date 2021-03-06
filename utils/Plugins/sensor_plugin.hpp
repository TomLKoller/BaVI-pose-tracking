/*
 * sensor_plugin.hpp
 *
 *  Created on: 25.07.2018
 *      Author: tomlucas
 */

#ifndef PLUGINS_SENSOR_PLUGIN_HPP_
#define PLUGINS_SENSOR_PLUGIN_HPP_

#include <vector>
#include <memory>
#include <mutex>

#include "ZaVI_Utils.hpp"
#include "threaded_object.hpp"
namespace zavi
::plugin {
	/**
	 * Base class for all sensor plugins
	 */
	class SensorPlugin : public ThreadedObject {
	public:
		double nextTime=-1;     //!< the timestamp of the next line
		double oldTime=0;//!< the old timestamp

		virtual ~SensorPlugin() {
			callbacks.clear();
		}
		SensorPlugin() {

		}

		/**
		 * adds a new callback to this imu
		 *
		 * This can be used to append different estimators to this plugin
		 * @param callback The callback to use it gets a reference to this pointer, to the calling estimator and the current program time
		 * @param estimator The Estimator which wants to add the callback
		 */
		template<typename plugin_type,typename caller_type>
		void addDataCallback(std::shared_ptr<caller_type> caller,
				std::function<void (plugin_type *, caller_type *, double ) > callback
		) {
			callbacks.push_back(std::make_shared<CallbackContainer<plugin_type,caller_type>>(callback, caller));
		}

		template<typename lambdaType,typename caller_type,typename plugin_type>
		void addDataCallback(std::shared_ptr<plugin_type>, std::shared_ptr<caller_type> caller,
				lambdaType callback
		) {
			callbacks.push_back(std::make_shared<LambdaContainer<lambdaType,caller_type,plugin_type>>(callback, caller));
		}

		template<typename callback_type>
		void addDataCallback(std::shared_ptr<void> caller,std::shared_ptr<callback_type> & callback) {
			callbacks.push_back(std::make_shared<SensorCallbackTypeErasure>(callback,caller));

		}



		/**
		 * Signaler, that new data is available
		 *
		 * Calls all registered callbacks. Should be called after each readed data point.
		 *
		 * Locks for thread safety if using multiple sensors
		 */
		void newDataAvailable(double time_diff) {
			mtx.lock();
			for (std::shared_ptr<SensorCallback> callback : callbacks) {
				callback->call(this, time_diff);
			}
			mtx.unlock();
		}

		/**
		 * The plugin reads the next state and stores it internally
		 * @return true if a next state was available, false otherwise
		 */
		virtual bool readNextState(){throw "Not implemented readNextState";}

		/**
		 *
		 * @return Whether this sensor_plugin can work batch_wise
		 */
		virtual bool can_batchwise() {return false;}

		/**
		 * Signaler, that batch is ready
		 *
		 * Calls all registered callbacks. Should be called after each readed data point.
		 *
		 * Locks for thread safety if using multiple sensors
		 */
		void batchReady() {
			batch_mtx.lock();
			double time = zavi::Timer::getNow();
			for (std::shared_ptr<SensorCallback> callback : batch_callbacks) {
				callback->call(this, time);
			}
			batch_mtx.unlock();
		}

		/**
		 * adds a new callback to this imu
		 *
		 * This can be used to append different estimators to this plugin
		 * @param callback The callback to use it gets a reference to this pointer, to the calling estimator and the current program time
		 * @param estimator The Estimator which wants to add the callback
		 */
		template<typename plugin_type,typename caller_type>
		void addBatchCallback(
				void (*callback)(plugin_type *, caller_type *, double ),
				std::shared_ptr<caller_type> caller) {
			assert(can_batchwise());     //Not a batch sensor
			batch_callbacks.push_back(std::make_shared<CallbackContainer<plugin_type,caller_type>>(callback, caller));
		}

		template<typename callback_type>
		void addBatchCallback(std::shared_ptr<callback_type> callback, std::shared_ptr<void> caller) {
			assert(can_batchwise());     //Not a batch sensor
			batch_callbacks.push_back(std::make_shared<SensorCallbackTypeErasure>(callback,caller));

		}

	private:
		static inline std::mutex mtx;     // mutex to block simultanous calls of newDataAvailable
		static inline std::mutex batch_mtx;// mutex to block simultanous calls of batchReady

		class SensorCallback {
		public:
			virtual ~SensorCallback() {}
			virtual void call(SensorPlugin * plugin, double time)const =0;
		};

		class SensorCallbackTypeErasure : public SensorCallback {
		public:

			struct TypeErasure {
				virtual ~TypeErasure() {}

				virtual void call(SensorPlugin * plugin, double time)const=0;

			};
			template<typename callback_type>
			struct CallbackWrapper: public TypeErasure {
				CallbackWrapper(std::shared_ptr<callback_type> callback, std::shared_ptr<void> caller):callback(callback),caller(caller) {}

				virtual void call(SensorPlugin * plugin, double time)const {
					callback->operator()(plugin, caller.get(), time);
				}
				std::shared_ptr<callback_type> callback;
				std::shared_ptr<void> caller;
			};

			virtual void call(SensorPlugin * plugin, double time)const {
				callback->call(plugin,time);
			}
			template<typename callback_type>
			SensorCallbackTypeErasure(std::shared_ptr<callback_type> callb, std::shared_ptr<void> &caller) :
			callback(new CallbackWrapper<callback_type>(callb,caller)) {
			}
			SensorCallbackTypeErasure(SensorCallbackTypeErasure &callb):
			callback(callb.callback) {
			}

		protected:
			std::shared_ptr<TypeErasure> callback;

		};

		/**
		 * Just a little container class for added callbacks
		 */
		template<typename plugin_type, typename caller_type>
		class CallbackContainer: public SensorCallback {
		private:
			std::function<void (plugin_type *, caller_type *, double ) > callback;     //< the callback function
			std::shared_ptr<caller_type> caller;//< the caller

		public:
			/**
			 * Construct a new container
			 * @param callback the callback to contain
			 * @param estimator the Estimator which uses the plugin
			 */
			CallbackContainer(
					std::function<void (plugin_type *, caller_type *, double ) > callback, std::shared_ptr<caller_type> caller) :
			callback(callback), caller(caller) {

			}
			/**
			 * call the contained callback
			 * @param plugin the plugin on which the callback is called
			 * @param time the timestamp when it is called
			 */
			virtual void call(SensorPlugin * plugin, double time_diff) const {
				callback((plugin_type *)plugin, caller.get(), time_diff);
			}
		};

/**
		 * Just a little container class for added callbacks
		 */
		template<typename lambdaType,typename caller_type,typename plugin_type>
		class LambdaContainer: public SensorCallback {
		private:
			lambdaType callback;     //< the callback function
			std::shared_ptr<caller_type> caller;//< the caller

		public:
			/**
			 * Construct a new container
			 * @param callback the callback to contain
			 * @param estimator the Estimator which uses the plugin
			 */
			LambdaContainer(
					lambdaType callback, std::shared_ptr<caller_type> caller) :
			callback(callback), caller(caller) {

			}
			/**
			 * call the contained callback
			 * @param plugin the plugin on which the callback is called
			 * @param time the timestamp when it is called
			 */
			virtual void call(SensorPlugin * plugin, double time_diff) const {
				callback((plugin_type*)plugin, caller.get(), time_diff);
			}
		};

		public:
		std::vector<std::shared_ptr<SensorCallback >> callbacks;     //< list of all appended callbacks
		std::vector<std::shared_ptr<SensorCallback >> batch_callbacks;//< list of all appended callbacks
	};

}

#endif /* PLUGINS_SENSOR_PLUGIN_HPP_ */
