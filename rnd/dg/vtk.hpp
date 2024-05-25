#ifndef __CFD24_GEOM_VTK_HPP__
#define __CFD24_GEOM_VTK_HPP__

#include "common.hpp"
#include <fstream>

struct VtkUtils{
	/// adds vtk caption
	static void append_header(std::string caption, std::ostream& fs);

	/// adds point list
	static void append_points(const std::vector<double>& points, std::ostream& fs);

	/**
	 * @brief adds cell data to saved vtk grid
	 * @param data      data vector
	 * @param data_cap  data caption as it should appear in vtk
	 * @param fname     file name of saved vtk grid
	 * @param ndata     length of data vector. If =0 uses data.size() as its length
	 */
	static void add_cell_data(const std::vector<double>& data,
	                          std::string data_cap,
	                          std::string fname,
	                          size_t ndata=0);

	/**
	 * @brief adds vertex data to the saved vtk file
	 * @param data      data vector
	 * @param data_cap  data caption as it should appear in vtk
	 * @param fname     saved vtk file name
	 * @param ndata     length of data vector. If =0 uses data.size() as its length
	 */
	static void add_point_data(const std::vector<double>& data,
	                           std::string data_cap,
	                           std::string fname,
	                           size_t ndata=0);

	/**
	 * @brief Writes data for increasing time points into the .vtk.series file
	 */
	class TimeSeriesWriter{
	public:
		/**
		 * @brief constructor
		 * @param stem file name stem
		 *
		 * This creates the series file named "<stem>.vtk.series".
		 * Vtk files with instant data will be saved into the "<stem>" directory.
		 * If this directory already exists it will be purged.
		 */
		TimeSeriesWriter(std::string stem);

		/**
		 * @brief adds new time point to the series file
		 * @param   tm time point value
		 * @returns file name in <stem> directory or empty string if the given time point is
		 *          not valid for save due to the time step condition
		 *
		 * This function only makes a record in the series file but does not create vtk file with
		 * instant data in the "<stem>" directory. The latter should be done manually using
		 * the returned filename.
		 */
		std::string add(double tm);

		/**
		 * @breif set saving time step
		 * @param tm_step step value
		 *
		 * Default step value is zero, that will make saving for each time point.
		 * If non zero value is set as a step then time save points will be calculated 
		 * only a single time point will be chosen for [N*tm_step-eps, (N+1)*tm_step] time period.
		 * TimeDependentWriter::add() calls with ignored time points will return empty strings.
		 */
		void set_time_step(double tm_step, double eps = 1e-6);
	private:
		const std::string _stem;
		const std::string _series_fn;
		std::string _fileslist;
		double _step = 0;
		double _step_eps = 0;
		int _last_saved_point = -1;
		void save_series() const;
		int get_time_point_index(double tm) const;
	};
};

#endif
