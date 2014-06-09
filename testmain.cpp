#include "sparsegrid/sparsegrid.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include "sparsegrid/filereading.hpp"

int main(){

	//------------------
	// Reading from file
	matrix data;
	data = read_csv_matrix("data.csv");
	//-----------------


	int ndim = 2;
	vector levels(ndim);
	levels.fill(5);
	bool full = false;

	RegularGrid a(ndim, levels, full);
	a.Initialize();

	a.EvaluateData(data);

	std::ofstream myfile;
	myfile.open("example.txt");
	a.ProjectData();
	myfile << a.EvalPointsGrid(21);
	myfile.close();
}