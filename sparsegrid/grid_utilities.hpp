//
// grid_utilities header file
//
// Collection of utilitiy functions for the grid structs. Defined a 'gridutil' 
// namespace to make it clear when these are being used.
//
// Author : Padarn Wilson
// Date : 19/06/14
//

#include "grid_typedefs.hpp"

namespace gridutils{

	//
	// PowInt: Takes the power a**b where both inputs and the output
	//         are known to be integers. 
	//
	int PowInt(int a,int b);

	//
	// RegularGridSize1D: Takes level and boundary condition to calculate
	//                    the one dimensional size of the grid.
	//
	int RegularGridSize1D(int level, bool boundary);

	//
	// RegularGridSize: Takes a vector of the levels in each dimension and
	//                  returns the total size of a regular grid with those
	//                  levels.
	//
	int RegularGridSize(vector level, bool boundary);

	//
	// RegularGridSizeVector: Takes a vector of the grid levels and returns a
	//                        new vector of the one dimensional sizes.
	//
	vector RegularGridSizeVector(vector level, bool boundary);

	//
	// StrideFromSize: Converts an vector of sizes in the corresponding strides
	//                 by shifting to the right and inserting 1.
	//
	vector StrideFromSize(vector sizes);

	//
	// IncreaseBit: Increases a 'bit vector'. This is really an internal 
	//				function for iterating over indices.
	// Input notes: bit - is modified in place
	//                    resets to zero if no increase possible
	//              dimsizes - is size in each dim to interatve over
	//              d - is a dimension to hold fixed (-1 for no hold)
	//              dvec - alternative version where vector of d is held fixed
	//
	void IncreaseBit(vector & bit, vector dimsizes, int d);
	void IncreaseBit(vector & bit, vector dimsizes, vector dvec);

	//
	// DristributeBit: Helper function to move bit indicies to respect
	//                 size maxes.
	//
	void DistributeBit(vector & bit, vector sizes);

	//
	// IncreaseBitMonotone: Increases a 'bit vector' corresponding to the 
	//                      grid indicies for d-dimensions, so that continued 
	//                      incrementing leads to monotone traversal of grid
	//
	void IncreaseBitMonotone(vector & bit, vector dimsizes);

	//
	// CornerStrides: Returns a vector corresponding to the strides in each
	//                dimension to index the 'left most corner index', which
	//                is the smallest grid index of those adjacent to x.
	// 
	vector CornerStrides(vector x, vector level, bool boundary);

	//
	// BoundaryBitAndFixed: Take the strides for a corner and creates a valid
	//						first bit and fixed vector. Modified inplace.
	//						returns 'false' if the point is outside of the
	//						domain.
	//
	bool BoundaryBitAndFixed(vector strides, vector corner, bool boundary, 
							 vector & bit, vector & fixed, int & nfixed);

	//
	// IndexX: Finds the x position of an index.
	//
	vector IndexX(vector index, vector level, bool boundary);

	//
	// RelativeX: Finds the relative x position between a point x and the x
	//            position of a given index. 
	//
	vector RelativeX(vector x, vector index, vector level, bool boundary);

	//
	// RelativeXScaled: Finds the relative x position between point x and the 
	//                  position an index scaled to [-1,1]
	//
	vector RelativeXScaled(vector x, vector index, vector level, bool boundary);

	// 
	// HatVal: Evaluates hat function at x - assumes [-1,1] domain.
	//
	double HatVal(vector x);

};