//
// Adaptive grid struct header file
//
// Author : Padarn Wilson
// Date : 01/07/14
//


#ifndef ADAPTIVEGRID_HPP
#define ADAPTIVEGRID_HPP

#include "Eigen/Dense"
#include "grid_typedefs.hpp"
#include <vector>
#include <queue>
#include <memory>
#include <tr1/unordered_set> 
#include "grid_operations.hpp"
#include "grid_utilities.hpp"
#include "combination_grid.hpp"
#include "lattice_tools.hpp"
#include "sub_grid.hpp"


namespace std {
namespace tr1 {
template<>
struct hash<Eigen::VectorXd>
{

    size_t make_hash(const int& v) const
    {
        return hash<int>()(v);
    }
    
    void hash_combine(std::size_t& h, const std::size_t& v) const
    {
        h ^= v + 0x9e3779b9 + (h << 6) + (h >> 2);
    }

    size_t operator()(const Eigen::VectorXd & v) const
    {
        size_t h=0;
        for(int i=0; i< v.size(); i++){
            hash_combine(h, make_hash((int) v(i)));
        }
        return h;
    }
    
};
}}

typedef std::pair<double, vector> vecpair;

class compare_pair
{
public:
  compare_pair(){};
  bool operator() (const vecpair& lhs, 
  				    const vecpair& rhs) const
  {
    return (lhs.first<rhs.first);
  }
};

typedef struct AdaptiveGrid
{
	matrix data_;
	int ndims_; 
	bool boundary_;
	std::tr1::unordered_set<vector> subgrid_levels_;
	std::priority_queue<vecpair,std::vector<vecpair>,compare_pair> queued_;

	AdaptiveGrid(int ndims, bool boundary);
	~AdaptiveGrid(){};
	void Initialize();
	void Refine();
	double CDFPriority(vector level);
	CombinationGrid * toCombination();

} AdaptiveGrid ;

#endif