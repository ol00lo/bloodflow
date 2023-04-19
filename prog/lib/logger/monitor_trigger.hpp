#ifndef BF_MONITOR_TRIGGER_HPP
#define BF_MONITOR_TRIGGER_HPP

#include "common.hpp"
#include <cmath>

namespace bf{

class IMonitorTrigger{
public:
	virtual ~IMonitorTrigger() = default;
	virtual bool check(size_t iter, double time) = 0;
};

class MonitorTrigger_Always: public IMonitorTrigger{
public:
	bool check(size_t, double) override{
		return true;
	}
};

class MonitorTrigger_Never: public IMonitorTrigger{
public:
	bool check(size_t, double) override{
		return false;
	}
};

class MonitorTrigger_EachNIter: public IMonitorTrigger{
public:
	MonitorTrigger_EachNIter(size_t iter_step): _iter_step(iter_step) {}
	bool check(size_t iter, double) override{
		return iter % _iter_step == 0;
	}
private:
	const size_t _iter_step;
};

class MonitorTrigger_TimePeriod: public IMonitorTrigger{
public:
	MonitorTrigger_TimePeriod(double time_step): _prev_save(0), _time_step(time_step) {}
	bool check(size_t, double time) override{
		size_t cur = (size_t)std::floor(time/_time_step);
		if (cur > _prev_save){
			_prev_save = cur;
			return true;
		} else {
			return false;
		}
	}
private:
	size_t _prev_save;
	const double _time_step;
};


}

#endif
