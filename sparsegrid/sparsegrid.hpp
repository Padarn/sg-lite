
/**
 * \file sparsegrid.hpp
 * \author Padarn Wilson
 * \date 05/06/14 
 * \brief Declaration for basic sparse grid structure.
 * \note1: 05/06/14 Starting with everything in this file to start with, trying to keep things light
 * 		   weight, and so will only decompose as needed. 
 * \note2: 18/06/14 - Have decided that I will store more in each 'regular' grid for ease of code 
 *         readability in the first instance. Can move things for performance later.
  */


#ifndef SPARSEGRID_HPP
#define SPARSEGRID_HPP

#include "Eigen/Dense"
#include <vector>
#include <memory>

/**
* Typedefs for clarity
*/
// TODO - Should add a way of compiling with fixed vector length
// NOTE - Even if this is possible the data vector will change
typedef Eigen::VectorXd vector;   
typedef Eigen::MatrixXd matrix;  

//
// RegualarGrid: Struct that defines a regular full anisotropic grid with 
//               hat basis functions.
//
typedef struct RegularGrid
{
	int ndims_;  // Number of dimensions.
	int size_;  // Total number of points in grid.
	vector levels_;  // Vector of levels in each dimension.
	vector strides_; // Vector giving stride in each dimension for 'data_'.
	vector sizes_; // Sizes in each dimension.
	vector data_; // Data vector - store single value at each point on grid.
	bool boundary_; // Flag to determine whether or not boundary exists.
	bool isInitialised_ = false; // Flag to check initialisation.

	// 
	// Constructor: Create a new Grid object. Data is not initialised at this point
	//              and so Initialise msut be called if it will be used.
	// Input
	// -----
	// ndims - dimensions of the grid
	// level - vector of level in each dimension
	// boundary - bool determining if boundary is included
	// 
	RegularGrid(int ndims, vector level, bool boundary);
	~RegularGrid(){};

	// 
	// Initialize: Allocates memory for the grid data array. Except for the case of 
	//             very small grids this will be the major memory use for the grid
	//             and so can remain uninitialised until needed.
	//
	void Initialize();

	//
	// EvaluateData: Evaluates a set of data points at each point in the grid and 
	//               stores the result averaged over all points.
	// Input
	// -----
	// data - eigen (const) matrix of data
	// or
	// datafile - string giving path to data file. Should be CSV format.
	//
	void EvaluateData(const matrix & data);  // matrix input
	void EvaluateData(std::string filename); // filename input

	//
	// ProjectionDensity: Treats the stored data vector of the grid as an evaluation
	//                    of point data over the grid (see EvaluateData) and then
	//                    calculates coefficients of density estimate based on L^2
	//                    projection.
	//
	void ProjectionDensity();

	// 
	// EvaluatatePoints: Takes a set of points and evaluates the full grid at those
	//                   points assuming that the stored data vector corresponds to
	//                   coefficients of the nodal basis
	// Input 
	// -----
	// data - matrix of data points
	//
	// Output
	// ------
	// vector of grid values at each of the input points
	//
	vector EvalPoints(matrix data);

	//
	// EvaluateGrid: Evaluates the grid over a regularly spaced grid with resolution
	//               as specified. Intended for plotting.
	// Input
	// -----
	// res - resolution of grid to use, e.g 10 gives 10 points in each dimension.
	//
	// Output
	// ------
	// vector of values over grid
	//
	vector EvalPointsGrid(int res);


	// EXPERIMENTAL ---------------------------------------
	void HeatSolve(double time);
	void CentralStepStart();
	// ----------------------------------------------------

} RegularGrid ;

// Typedef a pointer
typedef std::shared_ptr<RegularGrid> regulargridptr;


//
// CombinationGrid: Struct that defines a structure for a combination grid.
//                  The combination grid is made up of a set of different 
//                  regular grids and coefficients to combine them.
//
typedef struct CombinationGrid
{

	int ndims_;
	int boundary_;
	std::vector<vector> levels_;
	std::vector<double> coefs_; // float, int?
	std::vector<regulargridptr> grids_;



	CombinationGrid(int ndims, bool boundary);
	~CombinationGrid(){};
	// ---------------------
	//     BUILD GRIDS
	// ---------------------

	void Initialize();

	void SetupStandardCombinationGrid(int level);

	// ---------------------
	//   DO STUFF TO GRIDS
	// ---------------------

	// Evaluate data on grid
	void EvaluateData(const matrix & data);

	// L2 Projection of data (hat functions)
	void ProjectData();

	// ---------------------
	//     EVAL GRIDS
	// ---------------------

	vector EvalPoints(matrix data);

	vector EvalPointsGrid(int res);

	void HeatSolve(double time);

	void CentralStepStart();

} CombinationGrid;






#endif