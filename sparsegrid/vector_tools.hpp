//
// vector_tools header file
//
// Collection of tools for working with vector objects. Just quick helper 
// functions. Some of these may exist in Eigen, and I will try and simplify
// in the future
//
// Author : Padarn Wilson
// Date : 27/06/14
//

#include "grid_typedefs.hpp"
#include "combination_grid.hpp"

namespace vectortools{

	//
	// MinVector: Performs a componentwise minimum on two vectors.
	//
	vector MinVector(vector x, vector y);

	//
	// MaxVector: Performs a componentwise maximum on two vectors.
	//
	vector MaxVector(vector x, vector y);

	// 
	// VectorToIndex: Takes a vector and a max-size and converts it to
	// an index.
	//
	int VectorToIndex(vector in, vector maxlevel);

	//
	// IndexToVector: Does the opposite transformation to that of VectorToIndex
	//
	vector IndexToVector(int in, vector maxlevel);

};