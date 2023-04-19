#include "vtk_monitor.hpp"
#include <sstream>
#include <fstream>
#include <iomanip>

using namespace bf;

VtkMonitor::VtkMonitor(const Grid& grid, std::string fn, std::vector<std::string> fnames): _grid(grid), _fname(fn), _series_fname(fn+".vtk.series"), _data_names(fnames){
	_series_content.push_back("{");
	_series_content.push_back("  \"file-series-version\" : \"1.0\",");
	_series_content.push_back("  \"files\" : [");
}

VtkMonitor::VtkMonitor(const Grid& grid, std::string fn): VtkMonitor(grid, fn, {fn}){}

void VtkMonitor::save_series(std::string new_fn, double new_time){
	if (_series_content.size() > 3){
		_series_content.back().push_back(',');
	}
	std::ostringstream ofs;
	ofs << "    {\"name\": \"" << new_fn << "\", \"time\": " << new_time << "}";
	_series_content.push_back(ofs.str());

	std::ofstream fs(_series_fname);
	for (const std::string& s: _series_content){
		fs << s << std::endl;
	}
	fs << "  ]" << std::endl;
	fs << "} " << std::endl;
	fs.close();
}

void VtkMonitor::step(size_t iter, double time, std::vector<const std::vector<double>*> solution){
	// write data
	std::ostringstream oss;
	//oss << set_precision(3);
	oss << _fname << "_";
	oss << std::setfill('0') << std::setw(8) << iter;
	oss << "_";
	oss << std::setprecision(5) << time << ".vtk";
	_grid.save_data(_data_names, solution, oss.str());

	save_series(oss.str(), time);
}
