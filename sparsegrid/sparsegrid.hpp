
/**
 * \file sparsegrid.hpp
 * \author Padarn Wilson
 * \date 05/06/14 
 * \brief Declaration for basic sparse grid structure.
 * \notes Starting with everything in this file to start with, trying to keep things light
 * 		  weight, and so will only decompose as needed. 
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

/**
 * \struct regulargrid 
 * \brief Defined regulargrid element - isotropic regular grid.
 */
typedef struct RegularGrid
{

	int ndims_;
	int size_;
	vector levels_;
	vector data_;
	bool boundary_;

	// Basic Constructor and Destructor
	RegularGrid(int ndims, vector level, bool boundary);
	~RegularGrid(){};

	// Initialize the data vector
	void Initialize();

	// Evaluate data on grid
	void EvaluateData(const matrix & data);

	// L2 Projection of data (hat functions)
	void ProjectData();

	// Evaluate
	vector EvalPoints(matrix data);

	// Evaluate
	vector EvalPointsGrid(int res);


} RegularGrid ;

// Typedef a pointer
typedef std::shared_ptr<RegularGrid> regulargridptr;

/**
 * \struct regulargrid 
 * \brief Defined regulargrid element - isotropic regular grid.
 */
typedef struct SubGrid : RegularGrid
{

	// Basic Constructor and Destructor
	SubGrid(int ndims, vector level, bool boundary):
			 RegularGrid(ndims, level, boundary) {};
	~SubGrid(){};


} SubGrid ;

/**
 * \struct combinationgrid
 * \brief Collection of levels and coefficents. For now will have vector
 * of pointers to optional RegularGrids.
 */
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

} CombinationGrid;






#endif