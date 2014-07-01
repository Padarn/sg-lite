//
// CombinationGrid struct header file
//
// This file outlines the interface for the Combination grid. This implements a
// container that can be used to hold a combination of levels and perform
// operations across these 'components'. In the most basic set up, the 
// components are translated into RegularGrid structs with the classic
// combination grid coefficients. 
//
// Note: The storage is these levels is not sorted or indexed. Additional 
// structure is required ontop of this collection for that purpose.
//
// Author : Padarn Wilson
// Date : 26/06/14
//


#ifndef COMBINATIONGRID_HPP
#define COMBINATIONGRID_HPP

#include "Eigen/Dense"
#include "grid_typedefs.hpp"
#include "grid_operations.hpp"
#include "grid_utilities.hpp"
#include "file_reading.hpp"
#include "regular_grid.hpp"
#include <vector>
#include <memory>

//
// CombinationGrid: Struct that defines a structure for a combination grid.
//                  The combination grid is made up of a set of different 
//                  regular grids and coefficients to combine them.
//
typedef struct CombinationGrid
{

	int ndims_;
	bool boundary_;
	std::vector<vector> levels_;
	std::vector<double> coefs_; // float, int?
	std::vector<RegularGrid *> grids_;
	bool isInitialised_;


	// 
	// Constructor: Create a new empty combination container set up for an n 
	//              dimensional grid
	// 
	CombinationGrid(int ndims, bool boundary);
	~CombinationGrid();
	
	// 
	// Initialize: Allocates memory for each of the components. Currently just 
	//             regular grids.
	//
	void InitializeComponents();

	//
	// SetupStandardCombinationGrid: Sets up a standard combination grid at the
	//                               specified level.
	//
	void SetupStandardCombinationGrid(int level);


	/// ------------------------------ NOTE ----------------------------------
	/*
		The following functions are included for convienience at the time being
		I am not convinced these will be included in the end, instead perhaps
		preferring to manually do iteration over grid etc.
	*/
	void HatDataBin(const matrix & data);  // matrix input.
	void HatDataBin(std::string filename); // filename input.

	void HaarDataBin(const matrix & data); // matrix input.
	void HaarDataBin(std::string filename); // filename input.

	void CollectCDF();

	vector EvalPoints(matrix data);
	vector EvalPointsDerivative(matrix data);



} CombinationGrid;


#endif