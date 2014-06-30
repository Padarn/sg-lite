/*
*  Interface file for swig
*/

%module pysglite

%{
#define SWIG_FILE_WITH_INIT
#include <numpy/arrayobject.h>
#include "sparsegrid/regular_grid.hpp"
#include "sparsegrid/sub_grid.hpp"
#include "sparsegrid/file_reading.hpp"
#include "sparsegrid/combination_grid.hpp"
#include <Eigen/Core>
#include <Python.h>
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
%include "std_string.i"

%init
%{
  import_array();
%}

// Make swig aware of the typedefs
%typedef Eigen::VectorXd vector;   
%typedef Eigen::MatrixXd matrix; 

// Define fragment the does conversions
%fragment("EigenNumpyConversion", "header", fragment="NumPy_Fragments")
%{

Eigen::VectorXd * ConvertNumpyToEigenVector(PyObject* in)
{
 	int size = 0;
	 // check numpy array and get size
    PyObject *inarray = PyArray_FROM_OTF(in, NPY_DOUBLE, NPY_IN_ARRAY);
    if (inarray==NULL) std::cout << "ERROR: Conversion Failed." << std::endl;

    PyArrayObject *temp=NULL;
    if (PyArray_Check(inarray))
    {
        int newobj;
        temp = obj_to_array_allow_conversion(inarray, NPY_DOUBLE, &newobj);
        size = array_size(temp,0);
    }
  
    Eigen::VectorXd *out = new Eigen::VectorXd;
    (*out).resize(size);
    (*out).fill(0);
    double *  values = ((double *) PyArray_DATA(temp));
    for (long int i = 0; i != size; ++i){
        (*out)(i) = values[i];
    }
    return out;
}


Eigen::MatrixXd * ConvertNumpyToEigenMatrix(PyObject* in)
{
    int rows = 0;
    int cols = 0;
     // check numpy array and get size
    PyObject *inarray = PyArray_FROM_OTF(in, NPY_DOUBLE, NPY_IN_ARRAY);
    if (inarray==NULL) std::cout << "ERROR: Conversion Failed." << std::endl;

    // extract data
    PyArrayObject *temp=NULL;
    if (PyArray_Check(inarray))
    {
        int newobj;
        temp = obj_to_array_allow_conversion(inarray, NPY_DOUBLE, &newobj);
        rows = array_size(temp,0);
        cols = array_size(temp,1);
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

    void ConvertEigenToNumpyMatrix(PyObject** out, Eigen::MatrixXd * in)
    {
        npy_intp dims[2] = {in->rows(),in->cols()};
        *out = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
        double * data = static_cast<double*>(PyArray_DATA(*out));
        for (int i = 0; i != dims[0]; ++i){
            for (int j = 0; j != dims[1]; ++j){
                data[i*dims[1]+j] = (*in).coeff(i,j);
            }
        }
    };


%}


%typemap(in, fragment="EigenNumpyConversion") Eigen::VectorXd 
{
  $1 = *ConvertNumpyToEigenVector($input);
}

%typemap(in, fragment="EigenNumpyConversion") Eigen::VectorXd &
{
  $1 = ConvertNumpyToEigenVector($input);
}

%typemap(in, fragment="EigenNumpyConversion") Eigen::MatrixXd
{
  $1 = *ConvertNumpyToEigenMatrix($input);
}

%typemap(in, fragment="EigenNumpyConversion") Eigen::MatrixXd &
{ 
  $1 = ConvertNumpyToEigenMatrix($input);
}

%typemap(out, fragment="EigenNumpyConversion") Eigen::VectorXd
{
    ConvertEigenToNumpyVector(&$result, &$1);
}


%typemap(out, fragment="EigenNumpyConversion") Eigen::MatrixXd
{
    ConvertEigenToNumpyMatrix(&$result, &$1);
}

// Include the main header file with the typemaps defined
#include <vector>
%include "sparsegrid/regular_grid.hpp"
%include "sparsegrid/combination_grid.hpp"
%include "sparsegrid/sub_grid.hpp"


