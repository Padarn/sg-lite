
#ifndef FILE_FUNCTIONS_HPP
#define FILE_FUNCTIONS_HPP

#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>

typedef Eigen::MatrixXd matrix;  

double fromString( const std::string &s )
{
  double t;
  std::istringstream iss( s );
  iss >> t;
  return t;
}

int num_rows_in_file(const std::string& filename)
{
  int numLines = 0;
  std::string line;
  std::ifstream file;
  file.open( filename.c_str() );
  if ( file.is_open() )
    {
      while( !file.eof() )
	{
	  getline( file, line );
	  numLines++;
	}
      file.close();
      return numLines;
    }
  else
    {
      std::string msg = "num_rows_in_file() Unable to open file: ";
      msg += filename;
      throw( std::runtime_error( msg ) );
    }
}

int num_columns_in_file( const std::string& filename, char delimiter = ',' )
{
  std::string line;
  std::ifstream file;
  file.open( filename.c_str() );
  if ( file.is_open() )
    {
      getline(file, line);
      std::stringstream  lineStream( line );
      std::string        cell;
      int numCols = 0;
      while(std::getline( lineStream,cell, delimiter ) )
        {
	  numCols++;
        }
      file.close();
      return numCols;
    }
  else
    {
      std::string msg = "num_columns_in_file() Unable to open file: ";
      msg += filename;
      throw(std::runtime_error(msg));
    }
};


matrix read_csv_matrix( std::string filename,
	   char delimiter = ',' )
{
  int i,j;
  int numRows = num_rows_in_file( filename );
  int numCols = num_columns_in_file( filename );
  //Assumes the file is terminated with a "\n" as the last row
  matrix filecontents(numRows-1,numCols);
  filecontents.fill(0);
  
  std::ifstream file;
  file.open( filename.c_str()) ;
  if ( file.is_open() )
    {
      std::string line;
      i = 0;
      while( std::getline( file, line) )
	{
	  j = 0;
	  std::stringstream  lineStream( line );
	  std::string        cell;
	  while( std::getline( lineStream, cell, delimiter ) )
	    {
	      filecontents(i,j) = fromString( cell );
	      j++;
	    }
	  i++;
	}
      file.close();
    }
  else
    {
      std::string msg = "read( filename) Unable to open ";
      msg += "file " + filename;
      throw( std::runtime_error( msg ) );
    }
    return filecontents;

};

#endif