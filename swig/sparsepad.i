/*
*  Interface file for sparsepad
*  to compile:
*  swig -c++ -python -I../ sparsepad.i
*  g++ -std=c++11 -static-libstdc++ -fPIC -c sparsepad_wrap.cxx -o _sparsepad.o -L../sparsegrid -I/home/padarn/envs/statsmodelsdev/include/python2.7/ -I/home/padarn/envs/statsmodelsdev/lib/python2.7/site-packages/numpy/core/include/ -I/home/padarn/work/sparsepad -L/home/padarn/envs/statsmodelsdev/lib/python2.7/
*  g++ -shared ../sparsegrid/sparsegrid.o _sparsepad.o -o _sparsepad.so
*/

%module sparsepad

%{
#define SWIG_FILE_WITH_INIT
#include <Eigen/Core>
#include <Python.h>
#include <numpy/arrayobject.h>
#include "sparsegrid/regular_grid.hpp"
#include "sparsegrid/file_reading.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <memory>

typedef Eigen::VectorXd vector;   
typedef Eigen::MatrixXd matrix;  
%}

%include "numpy.i"
%include "eigen.i"
%include "std_string.i"

%typedef Eigen::VectorXd vector;   
%typedef Eigen::MatrixXd matrix; 

%fragment("PadarnTestFrag", "header", fragment="NumPy_Fragments")
%{
Eigen::VectorXd ConvertNumpyToEigenVector(PyObject* in)
{
 	int size = 0;
	 // check numpy array and get size
 	if (!is_array(in))
    {
      PyErr_SetString(PyExc_ValueError, "The given input is not known as a NumPy array or matrix.");
    }
 	// Check dimensions
    else if (array_numdims(in) != 1)
    {
      PyErr_SetString(PyExc_ValueError, "Wanted 1D Vector.");
    }
    else {
    	size = array_size(in,0);
    }

    // extract data
    int isNewObject = 0;
    PyArrayObject *temp=NULL;
    if (PyArray_Check(in))
    {
      temp = (PyArrayObject*)in;
    }

    Eigen::VectorXd out(size);
    out.fill(0);
    double *  values = ((double *) PyArray_DATA(temp));
    for (long int i = 0; i != size; ++i){
        out(i) = values[i];
    }
    return out;
}


Eigen::MatrixXd * ConvertNumpyToEigenMatrix(PyObject* in)
{
    int rows = 0;
    int cols = 0;
     // check numpy array and get size
    if (!is_array(in))
    {
      PyErr_SetString(PyExc_ValueError, "The given input is not known as a NumPy array or matrix.");
    }
    // Check dimensions
    else if (array_numdims(in) != 2)
    {
      PyErr_SetString(PyExc_ValueError, "Wanted 2D Matrix.");
    }
    else {
        rows = array_size(in,0);
        cols = array_size(in,1);
    }

    // extract data
    int isNewObject = 0;
    PyArrayObject *temp=NULL;
    if (PyArray_Check(in))
    {
      temp = (PyArrayObject*)in;
    }
    Eigen::MatrixXd *out = new Eigen::MatrixXd;
    (*out).resize(rows,cols);
    (*out).fill(0);
    double *  values = ((double *) PyArray_DATA(temp));
    for (long int i = 0; i != rows; ++i){
        for(long int j = 0; j != cols; ++j){
            (*out)(i,j) = values[i*rows+j];
        }
    }

    return out;

    }
   void ConvertEigenToNumpyVector(PyObject** out, Eigen::VectorXd * in)
   {
    npy_intp dims[1] = {in->size()};
    *out = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    double * data = static_cast<double*>(PyArray_DATA(*out));
    for (int i = 0; i != dims[0]; ++i){
        data[i] = (*in)(i);
      }
  };


%}


%typemap(in, fragment="PadarnTestFrag") Eigen::VectorXd
{
  $1 = ConvertNumpyToEigenVector($input);
}

%typemap(in, fragment="PadarnTestFrag", numinputs=0) Eigen::VectorXd & result
 ( Eigen::VectorXd TEMP)
{
  $1 = &TEMP;
};


%typemap(in, fragment="PadarnTestFrag") Eigen::MatrixXd
{
  $1 = *ConvertNumpyToEigenMatrix($input);
}

%typemap(in, fragment="PadarnTestFrag") Eigen::MatrixXd &
{ 
  $1 = ConvertNumpyToEigenMatrix($input);
}


%typemap(out, fragment="PadarnTestFrag") Eigen::VectorXd
{
    ConvertEigenToNumpyVector(&$result, &$1);
}



#include <vector>
%include "sparsegrid/regular_grid.hpp"


%init
%{
  import_array();
%}


