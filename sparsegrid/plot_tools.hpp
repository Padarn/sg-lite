//
// plot_tools header file
//
// Tools to simplify plotting.
//
// Author : Padarn Wilson
// Date : 24/06/14
//

#include "grid_typedefs.hpp"
#include "grid_utilities.hpp"

namespace plottools{

	//
	// MeshGrid: Creates a matrix of data points corresponding to a grid over
	//           [0,1]^d with 1D resolution specified by res.
	//
	matrix MeshGrid(int res, int d);

}