#include "sparsegrid/adaptive_grid.hpp"
#include "sparsegrid/plot_tools.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include "sparsegrid/file_reading.hpp"
#include <gtest/gtest.h>

int main(){

	//------------------
	// Reading from file
	matrix data;
	data = read_csv_matrix("data/data_heat.csv");
	//-----------------

	int ndim = 2;
	bool full = true;

	AdaptiveGrid a(ndim, full);
	a.data_ = data;
	a.Initialize();
	for(int i =0; i< 10; i++){
		a.Refine();
	}

	CombinationGrid * result = a.toCombination();
	for(int i = 0; i < result->levels_.size(); i++)
	{
		std::cout << "level " << result->levels_[i].transpose() << " coef " << 
		result->coefs_[i] << std::endl;
	}
	result->InitializeComponents();
	result->HaarDataBin(data);
	result->CollectCDF();

	matrix grid = plottools::MeshGrid(32,ndim);
	std::ofstream myfile;
	myfile.open("result.txt");
	myfile << result->EvalPointsDerivative(grid);
	myfile.close();

}