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
	data = read_csv_matrix("data_heat.csv");
	//-----------------


	int ndim = 2;
	vector levels(ndim);
	levels.fill(7);
	levels(1)=3;
	bool full = true;

	//RegularGrid a(ndim, levels, full);
	CombinationGrid a(ndim, full);
	a.SetupStandardCombinationGrid(8);
	a.Initialize();

	//a.CentralStepStart();
	a.EvaluateData(data);

	// std::ofstream myfile;
	// myfile.open("example1.txt");
	// a.ProjectData();
	// myfile << a.EvalPointsGrid(21);
	// myfile.close();

	// CombinationGrid b(ndim, full);
	// b.SetupStandardCombinationGrid(4);
	// b.Initialize();
	// b.EvaluateData(data);

	// std::ofstream myfile2;
	// myfile2.open("example2.txt");
	// b.ProjectData();
	// myfile2 << b.EvalPointsGrid(21);
	// myfile2.close();
	
	std::ofstream myfile2;
	myfile2.open("example_heat.txt");
	a.ProjectData();
	myfile2 << a.EvalPointsGrid(21);
	myfile2.close();

	std::ofstream myfile3;
	myfile3.open("example_heat3.txt");
	double t = 1.0;
	for (int i =0; i<2;i++){
		a.HeatSolve(t);
	}
	//a.ProjectData();
	myfile3 << a.EvalPointsGrid(31);
	myfile3.close();

}