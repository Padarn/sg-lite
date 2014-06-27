/*

TEST generalised_combination_test: Test will test the generalised combiantion
	 grid algorithm to make sure it works and produces the correct results;

*/

#include "sparsegrid/regular_grid.hpp"
#include "sparsegrid/combination_grid.hpp"
#include "sparsegrid/plot_tools.hpp"
#include "sparsegrid/lattice_tools.hpp"
#include "sparsegrid/file_reading.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>

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

int main(){


	int ndim = 3;
	bool full = true;
	int maxlevel = 4;

	CombinationGrid regular_combi(ndim, full);
	regular_combi.SetupStandardCombinationGrid(maxlevel);
	CombinationGrid adaptive_combi(ndim, full);

	// iterate through possible levels adding them one at a time
	vector onevec(ndim);// not sure I need this
	// interate bits and add for each level
	for(int i = 0; i <= maxlevel; i++)
	{	
		// calculate bits at this level
		int nlevel = nchoosek(ndim+i-1,i);
		onevec.fill(i);
		vector bit(ndim); bit.fill(0); bit(0)=i;
		// interate for nlevel and add
		for(int j = 0; j < nlevel; j++)
		{	
			latticetools::AddToCombination(&adaptive_combi,bit);
			gridutils::IncreaseBitMonotone(bit, onevec);
		}
	}

	// print results to check
	// TODO test using gtest
	std::cout << "normal combi" << std::endl;
	for(int i=0; i<regular_combi.levels_.size();i++)
	{
		vector level = regular_combi.levels_[i];
		double coef = regular_combi.coefs_[i];
		std::cout << "level " << level(0) << " " << level(1) << " " << level(2)
		          << " coef " << coef << std::endl;
	}
	std::cout << "adaptive combi" << std::endl;
	for(int i=0; i<adaptive_combi.levels_.size();i++)
	{
		vector level = adaptive_combi.levels_[i];
		double coef = adaptive_combi.coefs_[i];
		std::cout << "level " << level(0) << " " << level(1) << " " << level(2)
		          << " coef " << coef << std::endl;
	}

}