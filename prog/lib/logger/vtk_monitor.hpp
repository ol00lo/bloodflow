#ifndef BF_VTK_MONITOR_HPP
#define BF_VTK_MONITOR_HPP

#include "logger/imonitor.hpp"
#include "grid.hpp"

namespace bf{

class VtkMonitor: public IMonitor{
public:
	VtkMonitor(const Grid& grid, std::string fname);
	VtkMonitor(const Grid& grid, std::string fname, std::vector<std::string> data_names);
	void step(size_t iter, double time, std::vector<const std::vector<double>* > solution) override;
private:
	void save_series(std::string new_fn, double new_time);

	const Grid& _grid;
	std::string _fname;
	const std::string _series_fname;
	std::vector<std::string> _series_content;
	const std::vector<std::string> _data_names;
};

}
#endif
