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
	levels.fill(0);
	levels(0) = 3;
	bool full = true;

	RegularGrid a(ndim, levels, full);

	a.HatDataBin(data);
	a.CollectCDF();

	matrix grid = plottools::MeshGrid(31,ndim);
	std::ofstream myfile;
	myfile.open("result.txt");
	myfile << a.EvalPoints(grid);
	myfile.close();

}