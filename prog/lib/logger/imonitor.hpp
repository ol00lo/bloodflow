#ifndef BF_I_MONITOR_HPP
#define BF_I_MONITOR_HPP

#include "common.hpp"

namespace bf{

class IMonitor{
public:
	virtual ~IMonitor() = default;
	virtual void step(size_t iter, double time, std::vector<const std::vector<double>* > solution) = 0;

	void step(size_t iter, double time, const std::vector<double>& solution){
		step(iter, time, {&solution});
	}
};

}
#endif
