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
	bool full = true;

	RegularGrid a(ndim, levels, full);
	a.Initialize();

	a.EvaluateData(data);

	std::ofstream myfile;
	myfile.open("example1.txt");
	a.ProjectData();
	myfile << a.EvalPointsGrid(21);
	myfile.close();

	CombinationGrid b(ndim, full);
	b.SetupStandardCombinationGrid(4);
	b.Initialize();
	b.EvaluateData(data);

	std::ofstream myfile2;
	myfile2.open("example2.txt");
	b.ProjectData();
	myfile2 << b.EvalPointsGrid(21);
	myfile2.close();

}