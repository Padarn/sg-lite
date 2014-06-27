#include "lattice_tools.hpp"
#include <iostream>
#include <map>

namespace latticetools{

	CombinationGrid * MultiplyCombinations(CombinationGrid * grid1,
										   CombinationGrid * grid2)
	{
		int ndim = grid1->ndims_;		  // should check these two are the 
		bool boundary = grid1->boundary_; // same for grid1 and grid2
		CombinationGrid * result = new CombinationGrid(ndim, boundary);
		vector level_temp;

		for(int i = 0; i < grid1->levels_.size(); i++)
		{
			for(int j = 0; j < grid2->levels_.size(); j++)
			{
				// Multiplying levels results in minvector
				level_temp = vectortools::MinVector(grid1->levels_[i],
													grid2->levels_[j]);
				double coef = grid1->coefs_[i]*grid2->coefs_[j];

				// Add result to scheme
				result->levels_.push_back(level_temp);
				result->coefs_.push_back(coef);
			}
		}

		// There may be repetitions within the new Combination Grid.
		// Will need to get rid of these. Perhaps do when addding.

		return result;
	}
	
	void AddCombinations(CombinationGrid * grid1,
						 CombinationGrid * grid2)
	{
		// NOTE: This function may be slow. Once this is working should
		// rethink how the best way to do this is. Wait until building in
		// adaptive stuff, ass this may change how we want to do this.

		int ndim = grid1->ndims_;

		// Start by finding the maximum level in the two grids.
		vector maxvector(ndim);
		maxvector.fill(0);
		for(int i = 0; i < grid1->levels_.size(); i++)
		{
			maxvector = vectortools::MaxVector(maxvector, grid1->levels_[i]);
		}
		for(int i = 0; i < grid2->levels_.size(); i++)
		{
			maxvector = vectortools::MaxVector(maxvector, grid2->levels_[i]);
		}
		
		// Now using max level, create a set containing all of the levels
		// and the corresponding coefficients.
		std::map<int, double> coefs_set_;
		std::map<int, double>::iterator it;
		int setindex;
		for(int i = 0; i < grid1->levels_.size(); i++)
		{
			setindex = vectortools::VectorToIndex(grid1->levels_[i], maxvector);
			it = coefs_set_.find(setindex);
			if (it != coefs_set_.end())
			{
				(it->second)+=grid1->coefs_[i];
			} else {
				coefs_set_[setindex] = grid1->coefs_[i];
			}
		}
		for(int i = 0; i < grid2->levels_.size(); i++)
		{
			setindex = vectortools::VectorToIndex(grid2->levels_[i], maxvector);
			it = coefs_set_.find(setindex);
			if (it != coefs_set_.end())
			{
				(it->second)+=grid2->coefs_[i];
			} else {
				coefs_set_[setindex] = grid2->coefs_[i];
			}
		}
		
		// clear out CombinationGrids old values
		grid1->levels_.clear();
		grid1->coefs_.clear();

		// Add all levels in set to new levels with coefs
		for (it = coefs_set_.begin(); it != coefs_set_.end(); ++it)
		{
			if (it->second !=0) // dont add levels with zero coefficient
			{
	    		int index = it->first;
	    		vector level = vectortools::IndexToVector(index, maxvector);
	    		grid1->levels_.push_back(level);
	    		grid1->coefs_.push_back(it->second);
    		}
		}

	}

	void AddToCombination(CombinationGrid * grid, vector level)
	{
		
		if (grid->levels_.size()==0){

			grid->levels_.push_back(level);
			grid->coefs_.push_back(1);

		} else {
			// TODO: Seems finding the max level might be reasonably common
			// probably worth storing it somewhere.
			int ndims = grid->ndims_;
			bool boundary = grid->boundary_;

			vector maxvector = level;
			for(int i = 0; i < grid->levels_.size(); i++)
			{
				maxvector = vectortools::MaxVector(maxvector, grid->levels_[i]);
			}

			// build grid of changes - calculated by multiplying P_new by 
			// (1-P_old)
			// step 1) make P_new grid
			// step 2) multiply each existing level (1-P_old)
			CombinationGrid * changes = new CombinationGrid(ndims, boundary);
			CombinationGrid * null = new CombinationGrid(ndims, boundary);
			changes ->levels_.push_back(level);
			changes ->coefs_.push_back(1);
			for(int i = 0; i<grid->levels_.size(); i++)
			{
				CombinationGrid * temp = new CombinationGrid(ndims, boundary);
				temp->levels_.push_back(maxvector);
				temp->coefs_.push_back(1);
				temp->levels_.push_back(grid->levels_[i]);
				temp->coefs_.push_back(-1);
				CombinationGrid * temp2 = changes;
				changes = MultiplyCombinations(temp2, temp);
				AddCombinations(changes,null);
				delete temp;
				delete temp2;

			}
			// add to existing grid
			// step 3) add
			AddCombinations(grid, changes);

			// step 4) clean up
			delete changes;
			delete null;
		}
	}
}



