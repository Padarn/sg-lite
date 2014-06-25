/*

TEST cdf_haar_test: This test will run through the functionality of using the 
	 				cdf estimation on a full grid using the haar binning. The 
	 				result should be a linear interpolant of the histogram.

*/



#include "sparsegrid/regular_grid.hpp"
#include "sparsegrid/plot_tools.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include "sparsegrid/file_reading.hpp"

int main(){

	//------------------
	// Reading from file
	matrix data;
	data = read_csv_matrix("data/histsingle.csv");
	//-----------------

	int ndim = 2;
	vector levels(ndim);
	levels.fill(8);
	bool full = true;

	RegularGrid a(ndim, levels, full);

	a.HaarDataBin(data);
	a.CollectCDF();

	matrix grid = plottools::MeshGrid(31,ndim);
	std::ofstream myfile;
	myfile.open("result.txt");
	myfile << a.EvalPoints(grid);
	myfile.close();

}