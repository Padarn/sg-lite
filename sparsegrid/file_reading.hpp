
#ifndef FILE_FUNCTIONS_HPP
#define FILE_FUNCTIONS_HPP

#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>

typedef Eigen::MatrixXd matrix;  

double fromString( const std::string &s );

int num_rows_in_file(const std::string& filename);

int num_columns_in_file( const std::string& filename, char delimiter = ',' );


matrix read_csv_matrix( std::string filename,
	   char delimiter = ',' );

#endif