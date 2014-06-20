//
// grid_operations header file
//
// Collection of grid operations to perform on regular grids.
//
// Author : Padarn Wilson
// Date : 20/06/14
//

#include "grid_typedefs.hpp"
#include "grid_utilities.hpp"

namespace gridops {

	//
	// Projection1DSlice: Performes a projection along a one dimensional slice
	//                    using a tridiagonal solver.
	// Note: data is modified inplace
	//
	void Projection1DSlice(vector & data, int index_start, int index_stride,
						   int size, int level);

	//
	// Projection1D: Performed a projection along each slice in a single 
	//               dimension.
	// Note: data is modified in place.
	//
	void Projection1D(vector & data, vector strides, vector levels, int size, 
					  int dim);

	//
	// ProjectionND: Performs projection in each dimension sequentially.
	// Note: data is modified in place.
	//
	void ProjectionND(vector & data, vector strides, vector levels, int size);

}