/*

TEST generalised_combination_test: Test will test the generalised combiantion
	 grid algorithm to make sure it works and produces the correct results;

*/

#include "sparsegrid/regular_grid.hpp"
#include "sparsegrid/combination_grid.hpp"
#include "sparsegrid/plot_tools.hpp"
#include "sparsegrid/lattice_tools.hpp"
#include "sparsegrid/file_reading.hpp"
#include "sparsegrid/vector_tools.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

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


	// compare results
	// first collect into maps
	vector maxvector(ndim); maxvector.fill(maxlevel);
	std::map<int, double> normal_set;
	std::map<int, double> adaptive_set;
	std::map<int, double>::iterator it;
	std::map<int, double>::iterator it_a;

	std::cout << " level   normal   adaptive" <<std::endl;
	int setindex;
	for(int i = 0; i < regular_combi.levels_.size(); i++)
	{
		setindex = vectortools::VectorToIndex(regular_combi.levels_[i], maxvector);
		normal_set[setindex] = regular_combi.coefs_[i];
	}
	for(int i = 0; i < adaptive_combi.levels_.size(); i++)
	{
		setindex = vectortools::VectorToIndex(adaptive_combi.levels_[i], maxvector);
		it = adaptive_set.find(setindex);
		if (it != adaptive_set.end())
		{
			std::cout << "ERROR: index " << setindex << " is repeated in adaptive" << std::endl;
		} else {
			adaptive_set[setindex] = adaptive_combi.coefs_[i];
		}
	}
		
	// print out comparison for any in both
	for (it = normal_set.begin(); it != normal_set.end(); ++it)
	{
			int index = it->first;
    		vector level = vectortools::IndexToVector(index, maxvector);
    		std::cout << " " << level(0) << " " << level(1) << " " << level(2);
    		it_a = adaptive_set.find(index);
    		std::cout << "      " << it->second;
    		if (it_a != adaptive_set.end()){
    			std::cout << "       " << it_a->second << std::endl;
    		} else {
    			std::cout << "NA" << std::endl;
    		}
	}

	// look for any in the adaptive not in the regular combi
	for (it_a = adaptive_set.begin(); it_a != adaptive_set.end(); ++it_a)
	{
		int index = it_a->first;
    	vector level = vectortools::IndexToVector(index, maxvector);
	 	it = normal_set.find(index);
	 	if (it == normal_set.end()) std::cout << "ERROR: level " << level(0) << " " << 
	 	level(1) << " " << level(2) <<  " in adaptive but not normal." << std::endl;
	}
	
}