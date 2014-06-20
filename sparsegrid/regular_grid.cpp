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
	strides_ = gridutils::RegularGridSizeVector(levels_, boundary_);
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
		inDomain = gridutils::BoundaryBitAndFixed(strides_, corner, boundary_, 
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
				data_[bit.dot(index)] += gridutils::HatVal(relativeX); 
				gridutils::IncreaseBit(bit, strides_, fixed);
			}
		}
	}
}

vector RegularGrid::EvalPoints(matrix data)
{
	
	vector result(data.rows());
	result.fill(0);
	int ndata = data.rows();

	vector bit(ndims_);
	vector fixed(ndims_);
	vector relativeX(ndims_);
	vector index(ndims_);
	int nfixed;
	double value;
	bool inDomain;

	for(int i = 0; i < ndata; i++)
	{
		vector x = data.row(i);
		vector corner = gridutils::CornerStrides(x, levels_, boundary_);

		// get initial bit and fixed dimensions
		bit.fill(0); fixed.fill(0);
		inDomain = gridutils::BoundaryBitAndFixed(strides_, corner, boundary_, 
												  bit, fixed, nfixed);
		value = 0;
		// if x was in domain then it gives valid indicies so proceed
		if (inDomain)
		{
			
			int neighbours = gridutils::PowInt(2, ndims_ - nfixed);
			for(int j = 0; j < neighbours; j++)
			{
				
				index = bit+corner;
				relativeX = gridutils::RelativeXScaled(x, index,
													   levels_, boundary_);
				value += data_[bit.dot(index)]*gridutils::HatVal(relativeX);
				gridutils::IncreaseBit(bit, strides_, fixed);
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


