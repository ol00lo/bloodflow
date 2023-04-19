#ifndef BF_LOGGER_HPP
#define BF_LOGGER_HPP

#include "logger/imonitor.hpp"
#include "logger/monitor_trigger.hpp"

namespace bf{

class Logger{
public:
	void add_monitor(std::shared_ptr<IMonitor> monitor, std::shared_ptr<IMonitorTrigger> trigger);
	void step(size_t iter, double time, const std::vector<double>& solution);
	void step(size_t iter, double time, std::vector<const std::vector<double>*> solution);
private:
	std::vector<std::shared_ptr<IMonitor>> _monitors;
	std::vector<std::shared_ptr<IMonitorTrigger>> _triggers;
};

}
#endif
