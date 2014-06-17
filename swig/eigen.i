// -*- c++ -*-
 
%{
  #include <Eigen/Core>
%}

%include "numpy.i"


%fragment("Eigen_Fragments", "header",  fragment="NumPy_Fragments")
%{
  template <typename T> int NumPyType() {return -1;};
  
  template <class Derived>
  void ConvertFromNumpyToEigenMatrix(Eigen::MatrixBase<Derived>* out, PyObject* in)
  {
    int rows = 0;
    int cols = 0;
    // Check object type
    if (!is_array(in))
    {
      PyErr_SetString(PyExc_ValueError, "The given input is not known as a NumPy array or matrix.");
      return;
    }
    // Check data type
    else if (array_type(in) != NumPyType<typename Derived::Scalar>())
    {
      PyErr_SetString(PyExc_ValueError, "Type mismatch between NumPy and Eigen objects.");
      return;
    }
    // Check dimensions
    else if (array_numdims(in) > 2)
    {
      PyErr_SetString(PyExc_ValueError, "Eigen only support 1D or 2D array.");
      return;
    }
    else if (array_numdims(in) == 1)
    {
      rows = array_size(in,0);
      cols = 1;
      if ((Derived::RowsAtCompileTime != Eigen::Dynamic) && (Derived::RowsAtCompileTime != rows))
      {
        PyErr_SetString(PyExc_ValueError, "Row dimension mismatch between NumPy and Eigen objects (1D).");
        return;
      }
      else if ((Derived::ColsAtCompileTime != Eigen::Dynamic) && (Derived::ColsAtCompileTime != 1))
      {
        PyErr_SetString(PyExc_ValueError, "Column dimension mismatch between NumPy and Eigen objects (1D).");
        return;
      }
    }
    else if (array_numdims(in) == 2)
    {
      rows = array_size(in,0);
      cols = array_size(in,1);
      if ((Derived::RowsAtCompileTime != Eigen::Dynamic) && (Derived::RowsAtCompileTime != array_size(in,0)))
      {
        PyErr_SetString(PyExc_ValueError, "Row dimension mismatch between NumPy and Eigen objects (2D).");
        return;
      }
      else if ((Derived::ColsAtCompileTime != Eigen::Dynamic) && (Derived::ColsAtCompileTime != array_size(in,1)))
      {
        PyErr_SetString(PyExc_ValueError, "Column dimension mismatch between NumPy and Eigen objects (2D).");
        return;
      }
    }
    // Extract data
    int isNewObject = 0;
    PyArrayObject* temp = obj_to_array_contiguous_allow_conversion(in, array_type(in), &isNewObject);
    if (temp == NULL)
    {
      PyErr_SetString(PyExc_ValueError, "Impossible to convert the input into a Python array object.");
      return;
    }
    out->derived().setZero(rows, cols);
    typename Derived::Scalar* data = static_cast<typename Derived::Scalar*>(PyArray_DATA(temp));
    for (int i = 0; i != rows; ++i)
      for (int j = 0; j != cols; ++j)
        out->coeffRef(i,j) = data[i*cols+j];
  };
  
  template <class Derived>
  void ConvertFromEigenToNumPyMatrix(PyObject** out, Eigen::MatrixBase<Derived>* in)
  {
    npy_intp dims[2] = {in->rows(), in->cols()};
    *out = PyArray_SimpleNew(2, dims, NumPyType<typename Derived::Scalar>());
    typename Derived::Scalar* data = static_cast<typename Derived::Scalar*>(PyArray_DATA(*out));
    for (int i = 0; i != dims[0]; ++i)
      for (int j = 0; j != dims[1]; ++j)
        data[i*dims[1]+j] = in->coeff(i,j);
  };
  
  template<> int NumPyType<double>() {return PyArray_DOUBLE;};
%}

// ----------------------------------------------------------------------------
// Macro to create the typemap for Eigen classes
// ----------------------------------------------------------------------------
%define %eigen_typemaps(CLASS)

// In: (nothing: no constness)
%typemap(in, fragment="Eigen_Fragments") CLASS (CLASS temp)
{
  ConvertFromNumpyToEigenMatrix<CLASS>(&temp, $input);
  $1 = temp;
}

// In: Constructor
%typemap(in, fragment="Eigen_Fragments") CLASS (CLASS temp)
{
  ConvertFromNumpyToEigenMatrix<CLASS>(&temp, $input);
  $1 = temp;
}

// In: const&
%typemap(in, fragment="Eigen_Fragments") CLASS const& (CLASS temp)
{
  ConvertFromNumpyToEigenMatrix<CLASS>(&temp, $input);
  $1 = &temp;
}
// In: & (not yet implemented)
%typemap(in, fragment="Eigen_Fragments") CLASS &
{
  PyErr_SetString(PyExc_ValueError, "The input typemap for non-const reference is not yet implemented. Please report this problem to the developer.");
}
// In: const* (not yet implemented)
%typemap(in, fragment="Eigen_Fragments") CLASS const*
{
  PyErr_SetString(PyExc_ValueError, "The input typemap for const pointer is not yet implemented. Please report this problem to the developer.");
}
// In: * (not yet implemented)
%typemap(in, fragment="Eigen_Fragments") CLASS *
{
  PyErr_SetString(PyExc_ValueError, "The input typemap for non-const pointer is not yet implemented. Please report this problem to the developer.");
}

// Out: (nothing: no constness)
%typemap(out, fragment="Eigen_Fragments") CLASS
{
  ConvertFromEigenToNumPyMatrix<CLASS>(&$result, &$1);
}
// Out: const
%typemap(out, fragment="Eigen_Fragments") CLASS const
{
  ConvertFromEigenToNumPyMatrix<CLASS>(&$result, &$1);
}
// Out: const&
%typemap(out, fragment="Eigen_Fragments") CLASS const&
{
  ConvertFromEigenToNumPyMatrix<CLASS>(&$result, $1);
}
// Out: & (not yet implemented)
%typemap(out, fragment="Eigen_Fragments") CLASS &
{
  PyErr_SetString(PyExc_ValueError, "The output typemap for non-const reference is not yet implemented. Please report this problem to the developer.");
}
// Out: const* (not yet implemented)
%typemap(out, fragment="Eigen_Fragments") CLASS const*
{
  PyErr_SetString(PyExc_ValueError, "The output typemap for const pointer is not yet implemented. Please report this problem to the developer.");
}
// Out: * (not yet implemented)
%typemap(out, fragment="Eigen_Fragments") CLASS *
{
  PyErr_SetString(PyExc_ValueError, "The output typemap for non-const pointer is not yet implemented. Please report this problem to the developer.");
}

%enddef
