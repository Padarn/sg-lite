#include "combination_grid.hpp"
#include <iostream>

CombinationGrid::CombinationGrid(int ndims, bool boundary)
{
	ndims_ = ndims;
	boundary_ = boundary;
}


CombinationGrid::~CombinationGrid()
{
	if (isInitialised_)
	{
		for(int i = 0; i < grids_.size(); i++)
		{
			delete grids_[i];
		}
	}

}

void CombinationGrid::InitializeComponents()
{
	isInitialised_ = true;
	for(int i = 0; i < levels_.size(); i++)
	{
		RegularGrid * grid = new RegularGrid(ndims_, levels_[i], boundary_);
		grid->Initialize();
		grids_.push_back(grid);
	}
}

/// TODO: Move to a math utils class
int nchoosek(unsigned long n, unsigned long k) {
    unsigned long i;
    int b;
    if (0 == k || n == k) {
        return 1;
    }
    if (k > n) {
        return 0;
    }
    if (k > (n - k)) {
        k = n - k;
    }
    if (1 == k) {
        return n;
    }
    b = 1;
    for (i = 1; i <= k; ++i) {
        b *= (n - (k - i));
        if (b < 0) return -1; /* Overflow */
        b /= i;
    }
    return b;
}

void CombinationGrid::SetupStandardCombinationGrid(int level)
{ 
	/// NOTE :: CURRENT ASSUMING BOUNDARY

	vector onevec(ndims_);// not sure I need this
	// interate bits and add for each level
	for(int i = level-ndims_+1; i <= level; i++)
	{	
		// calculate bits at this level
		int nlevel = nchoosek(ndims_+i-1,i);
		// calculate coefficient at this level
		double coef = pow(-1,level-i)*nchoosek(ndims_-1,level-i);
		// get starting level
		onevec.fill(i);
		vector bit(ndims_); bit.fill(0); bit(0)=i;
		// interate for nlevel and add
		for(int j = 0; j < nlevel; j++)
		{
			levels_.push_back(bit);
			coefs_.push_back(coef);
			gridutils::IncreaseBitMonotone(bit, onevec);
		}
	}
}


// ------------------------- TEMP - SEE HEADER -----------------------------

void CombinationGrid::HatDataBin(std::string datafile)
{
	const matrix data = read_csv_matrix(datafile);
	HatDataBin(data);
}

void CombinationGrid::HatDataBin(const matrix & data)
{
	if (!isInitialised_) InitializeComponents();
	for(int i = 0; i < levels_.size(); i++)
	{
		RegularGrid * grid = grids_[i];
		if (!grid->isInitialised_) grid->Initialize();
		grid->HatDataBin(data);
	}
}  

void CombinationGrid::HaarDataBin(std::string datafile)
{
	const matrix data = read_csv_matrix(datafile);
	HaarDataBin(data);
}

void CombinationGrid::HaarDataBin(const matrix & data)
{
	if (!isInitialised_) InitializeComponents();
	for(int i = 0; i < levels_.size(); i++)
	{
		RegularGrid * grid = grids_[i];
		if (!grid->isInitialised_) grid->Initialize();
		grid->HaarDataBin(data);
	}
}  

void CombinationGrid::CollectCDF()
{
	if (!isInitialised_) InitializeComponents();
	for(int i = 0; i < levels_.size(); i++)
	{
		RegularGrid * grid = grids_[i];
		if (!grid->isInitialised_) grid->Initialize();
		grid->CollectCDF();
	}
}

vector CombinationGrid::EvalPoints(matrix data)
{
	if (!isInitialised_) InitializeComponents();
	vector result(data.rows());
	for(int i = 0; i < levels_.size(); i++)
	{
		RegularGrid * grid = grids_[i];
		if (!grid->isInitialised_) grid->Initialize();
		result += grid->EvalPoints(data) * coefs_[i];
	}
	return result;
}  