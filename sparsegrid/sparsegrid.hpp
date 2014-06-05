
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

/**
* Typedefs for clarity
*/
// TODO - Should add a way of compiling with fixed vector length
typedef Eigen::VectorXd vector;   


/**
 * \class subgrid 
 * \brief Defined subgrid element - either an isotropic regular grid
 * or a subspace of a full grid.
 */
struct subgrid
{

	int ndims;
	bool full;
	vector level;
	



}subgrid;


#endif