#include "adaptive_grid.hpp"
#include <iostream>

AdaptiveGrid::AdaptiveGrid(int ndims, bool boundary)
{
	ndims_ = ndims;
	boundary_ = boundary;
}

void AdaptiveGrid::Initialize()
{
	vector initial(ndims_);
	if (boundary_) initial.fill(0); 
	else initial.fill(1);
	queued_.push(vecpair(1.0,initial));
}

double AdaptiveGrid::CDFPriority(vector level){
	SubGrid test(ndims_, level);
	test.HaarDataBin(data_);
	return test.data_.lpNorm<Eigen::Infinity>();
}

//
// Basic Description - Add a new level and add levels to priority queue
//
void AdaptiveGrid::Refine()
{
	vector top_priority = queued_.top().second;
	queued_.pop();
	subgrid_levels_.insert(top_priority);
	int d = ndims_;
	vector bit(d);
	for(int i = 0; i<d; i++)
	{
		bit.fill(0);
		bit(i)=1;
		vector newvec = bit+top_priority;
		if (subgrid_levels_.find(top_priority+bit)==subgrid_levels_.end()){
			double priority = CDFPriority(newvec);
			queued_.push(vecpair(priority,newvec));
		}
	}
}

CombinationGrid * AdaptiveGrid::toCombination()
{
	CombinationGrid * result = new CombinationGrid(ndims_, true);
	for (auto it = subgrid_levels_.begin(); 
		it != subgrid_levels_.end(); ++it) {
    	latticetools::AddToCombination(result,*it);
		}
	return result;
}

