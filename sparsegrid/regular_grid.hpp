//
// RegularGrid struct header file
//
// This file outlines the interface for the RegularGrid struct which implements a
// regular anisotropic grid structure based on hat functions. Methods are are
// included for density estimation.
//
// Author : Padarn Wilson
// Date : 19/06/14
//


#ifndef REGULARGRID_HPP
#define REGULARGRID_HPP

#include "Eigen/Dense"
#include "grid_typedefs.hpp"
#include <vector>
#include <memory>
#include "grid_operations.hpp"
#include "grid_utilities.hpp"

//
// RegualarGrid: Struct that defines a regular full anisotropic grid with 
//               hat basis functions.
//
typedef struct RegularGrid
{
	int ndims_; 
	int size_;
	int ndata_;
	vector levels_;  
	vector strides_; 
	vector data_;
	bool boundary_;
	bool isInitialised_;

	// 
	// Constructor: Create a new Grid object. Data is not initialised at this
	//              point and so Initialise msut be called if it will be used.
	// 
	RegularGrid(int ndims, vector level, bool boundary);
	~RegularGrid(){};

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
	// ProjectionDensity: Treats the stored data vector of the grid as an 
	//                    evaluation of point data over the grid and then
	//                    calculates coefficients of density estimate based on 
	//                    L^2 projection.
	//
	void ProjectionDensity();

	// 
	// EvaluatatePoints: Takes a set of points and evaluates the full grid at 
	//                   those points assuming that the stored data vector 
	//                   corresponds to coefficients of the nodal basis
	//
	vector EvalPoints(matrix data);


} RegularGrid ;


#endif