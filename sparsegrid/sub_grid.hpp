//
// SubGrid struct header file
//
// This file outlines the interface for the SubGrid struct which implements 
// hierarchical subgrid.
//
// Author : Padarn Wilson
// Date : 30/06/14
//


#ifndef SUBGRID_HPP
#define SUBGRID_HPP

#include "Eigen/Dense"
#include "file_reading.hpp"
#include "grid_typedefs.hpp"
#include <vector>
#include <memory>
#include "grid_operations.hpp"
#include "grid_utilities.hpp"

//
// SubGrid: Struct that defines a regular hierarchical subgrid.
//
typedef struct SubGrid
{
	int ndims_; 
	int size_;
	int ndata_;
	vector levels_;
	vector sizes_; 
	vector strides_; 
	vector data_;
	bool isInitialised_;

	// 
	// Constructor: Create a new Grid object. Data is not initialised at this
	//              point and so Initialise msut be called if it will be used.
	// 
	SubGrid(int ndims, vector level);
	~SubGrid(){};

	// 
	// Initialize: Allocates memory for the grid data array. Except for the case 
	//             of very small grids this will be the major memory use for the
	//             grid and so can remain uninitialised until needed.
	//
	void Initialize();

	//
	// HatDataBin: Evaluates hat basis function at each grid point for set of 
	//             input data. Does not divide by ndata, but stores it instead.
	// Input Note: Input can either be the matrix of data or a filename.
	// Note: If data_ is not zero, this just adds to it.
	//
	void HatDataBin(const matrix & data);  // matrix input
	void HatDataBin(std::string filename); // filename input


	//
	// HaarDataBin: Evaluates haar basis function at each grid point for set of 
	//             input data. Does not divide by ndata, but stores it instead.
	// Input Note: Input can either be the matrix of data or a filename.
	// Note: If data_ is not zero, this just adds to it.
	//
	void HaarDataBin(const matrix & data); // matrix input
	void HaarDataBin(std::string filename); // filename input

	// 
	// EvaluatatePoints: Takes a set of points and evaluates the full grid at 
	//                   those points assuming that the stored data vector 
	//                   corresponds to coefficients of the nodal basis
	//
	vector EvalPoints(matrix & data);


} SubGrid;

#endif