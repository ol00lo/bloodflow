#include "logger.hpp"

using namespace bf;

void Logger::add_monitor(std::shared_ptr<IMonitor> monitor, std::shared_ptr<IMonitorTrigger> trigger){
	_monitors.push_back(monitor);
	_triggers.push_back(trigger);
}

void Logger::step(size_t iter, double time, std::vector<const std::vector<double>*> solution) {
	for (size_t i=0; i<_monitors.size(); ++i){
		if(_triggers[i]->check(iter, time)){
			_monitors[i]->step(iter, time, solution);
		}
	}
}

void Logger::step(size_t iter, double time, const std::vector<double>& solution) {
	step(iter, time, {&solution});
}
