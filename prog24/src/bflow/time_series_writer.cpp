#include "time_series_writer.hpp"
#include "ghc/filesystem.hpp"
#include <cmath>
#include <iomanip>
#include <sstream>

namespace fs = ghc::filesystem;
void create_directory(std::string path, bool purge)
{
    if (fs::is_directory(path))
    {
        if (purge)
        {
            fs::remove_all(path);
        }
        else
        {
            return;
        }
    }
    fs::create_directory(path);
}

TimeSeriesWriter::TimeSeriesWriter(std::string stem) : _stem(stem), _series_fn(stem + ".vtk.series")
{
    create_directory(stem, true);

    std::ofstream ofs(_series_fn);
    if (!ofs)
        throw std::runtime_error("Failed to open " + _series_fn + " for writing");

    ofs << "{" << std::endl;
    ofs << "  \"file-series-version\" : \"1.0\"," << std::endl;
    ofs << "  \"files\" : [" << std::endl;
    ofs << "  ]" << std::endl;
    ofs << "}" << std::endl;
}

std::string TimeSeriesWriter::add(double tm)
{
    if (_step > 0)
    {
        int itime_point = get_time_point_index(tm);
        if (itime_point > _last_saved_point)
        {
            _last_saved_point = itime_point;
        }
        else
        {
            return "";
        }
    }

    std::ostringstream fn;
    fn << std::setfill('0') << std::setw(8) << std::fixed << std::setprecision(4) << tm << ".vtk";
    std::string ret = _stem + '/' + fn.str();

    std::ostringstream oss;
    if (!_fileslist.empty())
    {
        oss << "," << std::endl;
    }
    oss << "    {\"name\": \"" << ret << "\", \"time\": " << tm << "}";
    _fileslist += oss.str();

    save_series();

    return ret;
}

void TimeSeriesWriter::save_series() const
{
    std::fstream ofs(_series_fn);
    if (!ofs)
        throw std::runtime_error("Failed to open " + _series_fn + " for writing");

    ofs << "{" << std::endl;
    ofs << "  \"file-series-version\" : \"1.0\"," << std::endl;
    ofs << "  \"files\" : [" << std::endl;
    ofs << _fileslist << std::endl;
    ofs << "  ]" << std::endl;
    ofs << "}" << std::endl;
}

int TimeSeriesWriter::get_time_point_index(double tm) const
{
    int index = int(std::floor(tm / _step));
    double cur_point = index * _step;
    double next_point = cur_point + _step;
    if (next_point - tm < _step_eps)
    {
        return index + 1;
    }
    else
    {
        return index;
    }
}

void TimeSeriesWriter::set_time_step(double tm_step, double eps)
{
    _step = tm_step;
    _step_eps = eps;
}
