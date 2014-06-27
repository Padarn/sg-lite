//
// lattice_tools header file
//
// Collection of tools for working with 'lattices'. The main motivation being
// the generalized sparse grid algorithm for building the combination method
// for adaptive sparse grids out of the subspaces included.
//
// NOTES: For the time being these methods are not implemented with efficiency
// in mind. I doubt this will be a bottlekneck, and for the time being clarity
// of algorithms and flexibility seem more important here.
//
// Author : Padarn Wilson
// Date : 26/06/14
//

#include "grid_typedefs.hpp"
#include "combination_grid.hpp"
#include "vector_tools.hpp"

namespace latticetools{

	//
	// MultiplyCombinations: This function takes two combination grid set ups
	//         consisting of levels and coefficents ant returns a new a new
	//         combination grid corresponding to the 'product' of the two rules.
	//         Product is to be understood in terms of a product of projections.
	// NOTE: Result may have duplicate elements.
	//         
	CombinationGrid * MultiplyCombinations(CombinationGrid * grid1,
										   CombinationGrid * grid2);

	//
	// AddCombinations: This function takes two combination grid set ups and 
	//                   calculates their sum - result stored in grid1
	// NOTE: Input may have duplicates, but result removes them.
	//
	void AddCombinations(CombinationGrid * grid1, CombinationGrid * grid2);

	//
	// AddToCombination: This function takes a combination grid (need not be
	//         the standard combination grid) and adds a new level to it. This
	//         produces a new set of levels and coefficients.
	//
	void AddToCombination(CombinationGrid * grid, vector level);

};