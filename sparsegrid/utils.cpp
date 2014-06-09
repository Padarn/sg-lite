#include "sparsegrid.hpp"
#include <iostream>

typedef Eigen::VectorXd vector;   
typedef Eigen::MatrixXd matrix;

void printvec(vector vec, std::string name)
{
	std::cout << name << std::endl;
	std::cout << vec.transpose() << std::endl;
}
