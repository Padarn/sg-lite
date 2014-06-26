/*

TEST combination_cdf_haar_test: This test will run through the functionality of 
					using the cdf estimation on a combination grid using the 
					haar binning. The result should be a linear interpolant of 
					the histogram.

*/

#include "sparsegrid/regular_grid.hpp"
#include "sparsegrid/combination_grid.hpp"
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
	bool full = true;

	CombinationGrid a(ndim, full);
	a.SetupStandardCombinationGrid(4);

	a.HaarDataBin(data);
	a.CollectCDF();

	matrix grid = plottools::MeshGrid(31,ndim);
	std::ofstream myfile;
	myfile.open("result.txt");
	myfile << a.EvalPoints(grid);
	myfile.close();

}