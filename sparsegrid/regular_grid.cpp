/* NOTES
*
* Useful article on indexing: http://nadeausoftware.com/articles/2012/06/c_c_tip_how_loop_through_multi_dimensional_arrays_quickly
*/

#include "regular_grid.hpp"
#include "file_reading.hpp"
#include "grid_utilities.hpp"
#include <iostream>

RegularGrid::RegularGrid(int ndims, vector levels, bool boundary)
{
	ndims_ = ndims;
	levels_ = levels;
	boundary_ = boundary;
	ndata_ = 0;
	size_ = gridutils::RegularGridSize(levels_, boundary_);
	sizes_ = gridutils::RegularGridSizeVector(levels_, boundary_);
	strides_ = gridutils::StrideFromSize(sizes_);
	isInitialised_ = false;
}

void RegularGrid::Initialize()
{
	data_ = vector(size_);
	data_.fill(0);
	isInitialised_ = true;
}

void RegularGrid::HatDataBin(std::string datafile)
{
	const matrix data = read_csv_matrix(datafile);
	HatDataBin(data);
}

void RegularGrid::HatDataBin(const matrix & data)
{
	if (!isInitialised_) Initialize();
	int ndata = data.rows();
	ndata_ += ndata;

	vector bit(ndims_);
	vector fixed(ndims_);
	vector relativeX(ndims_);
	vector index(ndims_);
	int nfixed;
	bool inDomain;

	for(int i = 0; i < ndata; i++)
	{
		vector x = data.row(i);
		vector corner = gridutils::CornerStrides(x, levels_, boundary_);

		// get initial bit and fixed dimensions
		bit.fill(0); fixed.fill(0);
		inDomain = gridutils::BoundaryBitAndFixed(sizes_, corner, boundary_, 
												  bit, fixed, nfixed);
		// if x was in domain then it gives valid indicies so proceed
		if (inDomain)
		{
			int neighbours = gridutils::PowInt(2, ndims_ - nfixed);
			for(int j = 0; j < neighbours; j++)
			{
				
				index = bit+corner;
				relativeX = gridutils::RelativeXScaled(x, index,
													   levels_, boundary_);
				data_[strides_.dot(index)] += gridutils::HatVal(relativeX); 
				gridutils::IncreaseBit(bit, sizes_, fixed);
			}
		}
	}
}

void RegularGrid::HaarDataBin(std::string datafile)
{
	const matrix data = read_csv_matrix(datafile);
	HatDataBin(data);
}

void RegularGrid::HaarDataBin(const matrix & data)
{
	if (!isInitialised_) Initialize();
	
	// TODO: Currently only handling situation for grid with boundary.
	if (!boundary_)
	{
		std::cout << "HaarDataBin Requires Boundary: No Action" << std::endl;
	}

	int ndata = data.rows();
	ndata_ += ndata;

	vector bit(ndims_);
	vector fixed(ndims_);
	vector index(ndims_);
	int nfixed;
	bool inDomain;

	for(int i = 0; i < ndata; i++)
	{
		vector x = data.row(i);
		vector corner = gridutils::CornerStrides(x, levels_, boundary_);

		// get initial bit and fixed dimensions
		bit.fill(0); fixed.fill(0);

		// NOTE: n should only be fixed if it is on the right boundary. So 
		// turn all other bits on.
		inDomain = gridutils::BoundaryBitAndFixed(sizes_, corner, boundary_, 
												  bit, fixed, nfixed);
		bit = 1 - fixed.array();
		// if x was in domain then it gives valid indicies so proceed
		if (inDomain)
		{		
			index = bit+corner;
			data_[strides_.dot(index)] += 1.0; 
		}
	}
}

vector RegularGrid::EvalPoints(matrix & data)
{
	// If uninitialized, assume constant.
	if (!isInitialised_)
	{
		Initialize();
		data_.fill(1.0);
	}
	vector result(data.rows());
	result.fill(0);
	int ndata = data.rows();

	vector bit(ndims_);
	vector fixed(ndims_);
	vector relativeX(ndims_);
	vector index(ndims_);
	vector onevec(ndims_); onevec.fill(1.0);
	int nfixed;
	double value;
	bool inDomain;

	for(int i = 0; i < ndata; i++)
	{
		vector x = data.row(i);
		vector corner = gridutils::CornerStrides(x, levels_, boundary_);

		// get initial bit and fixed dimensions
		bit.fill(0); fixed.fill(0); nfixed = 0;
		inDomain = gridutils::BoundaryBitAndFixed(sizes_, corner, boundary_, 
												  bit, fixed, nfixed);
		value = 0;
		// if x was in domain then it gives valid indicies so proceed
		if (inDomain)
		{
			
			int neighbours = gridutils::PowInt(2, ndims_ - nfixed);
			for(int j = 0; j < neighbours; j++)
			{
				
				index = bit+corner;
				vector indexx = gridutils::IndexX(index, levels_, boundary_);
				relativeX = gridutils::RelativeXScaled(x, index,
													   levels_, boundary_);
				value += data_[strides_.dot(index)]*gridutils::HatVal(relativeX);
				gridutils::IncreaseBit(bit, onevec, fixed);

			}
		}
		result(i)=value;
	}
	return result;

}

// Do L2 projection of stored data
void RegularGrid::ProjectionDensity()
{	
	gridops::ProjectionND(data_, strides_, levels_, size_);
}


