#include "sub_grid.hpp"
#include <iostream>

SubGrid::SubGrid(int ndims, vector level)
{
	ndims_ = ndims;
	levels_ = level;
	ndata_ = 0;
	size_ = gridutils::SubGridSize(levels_);
	sizes_ = gridutils::SubGridSizeVector(levels_);
	strides_ = gridutils::StrideFromSize(sizes_);
	isInitialised_ = false;
}


void SubGrid::Initialize()
{
	data_ = vector(size_);
	data_.fill(0);
	ndata_ = 0;
	isInitialised_ = true;
}

void SubGrid::HatDataBin(std::string filename)
{
	const matrix data = read_csv_matrix(filename);
	HatDataBin(data);
}


void SubGrid::HatDataBin(const matrix & data)
{
	if (!isInitialised_) Initialize();
	int ndata = data.rows();
	ndata_ += ndata;
	
	// count and locate boundary indices
	int nboundary = 0;
	vector indexboundary(ndims_); // end of vector will be unused
	for(int i = 0; i < ndims_; i++)
	{
		if (levels_(i)==0)
		{
			indexboundary(i) = 1; // dfixed vector
			nboundary++;
		}
	}
	vector onevec(ndims_); onevec.fill(1.0);

	// go through data points
	for(int i = 0; i < ndata; i++)
	{
		vector x = data.row(i);
		vector index = gridutils::SubGridIndex(x, levels_);
		// if no boundary
		if(nboundary == 0)
		{
			int bindex = index.dot(strides_);
			vector relativeX = gridutils::RelativeXScaledSubGrid(x, index, 
																 levels_);
			data_[bindex] += gridutils::HatVal(relativeX);
		}
		// else need to visit each boundary basis function
		else
		{
			vector bit(ndims_); bit.fill(0);
			for(int j =0; j < gridutils::PowInt(2,nboundary); j++)
			{
				int bindex = strides_.dot(index+bit);
				vector relativeX = gridutils::RelativeXScaledSubGrid(x, index, 
																	 levels_);
				data_[bindex] += gridutils::HatVal(relativeX);
				gridutils::IncreaseBit(bit, onevec, indexboundary);
			}
		}
	}

} 

void SubGrid::HaarDataBin(std::string filename)
{
	const matrix data = read_csv_matrix(filename);
	HaarDataBin(data);
}

void SubGrid::HaarDataBin(const matrix & data)
{
	if (!isInitialised_) Initialize();
	int ndata = data.rows();
	ndata_ += ndata;
	
	// count and locate boundary indices
	int nboundary = 0;
	vector indexboundary(ndims_); // end of vector will be unused
	for(int i = 0; i < ndims_; i++)
	{
		if (levels_(i)==0)
		{
			indexboundary(i) = 1; // dfixed vector
			nboundary++;
		}
	}
	vector onevec(ndims_); onevec.fill(1.0);

	// go through data points
	for(int i = 0; i < ndata; i++)
	{
		vector x = data.row(i);
		vector index = gridutils::SubGridIndex(x, levels_);
		// if no boundary
		if(nboundary == 0)
		{
			int bindex = index.dot(strides_);
			vector relativeX = gridutils::RelativeXScaledSubGrid(x, index, 
																 levels_);
			data_[bindex] += gridutils::HaarVal(relativeX);
		}
		// else need to visit each boundary basis function
		else
		{
			vector bit(ndims_); bit.fill(0);
			for(int j =0; j < gridutils::PowInt(2,nboundary); j++)
			{
				int bindex = strides_.dot(index+bit);
				vector relativeX = gridutils::RelativeXScaledSubGrid(x, index, 
																	 levels_);
				data_[bindex] += gridutils::HaarVal(relativeX);
				gridutils::IncreaseBit(bit, onevec, indexboundary);
			}
		}
	}

} 


vector SubGrid::EvalPoints(matrix & data)
{
	if (!isInitialised_) Initialize();
	int ndata = data.rows();
	ndata_ += ndata;
	
	// count and locate boundary indices
	int nboundary = 0;
	vector indexboundary(ndims_); // end of vector will be unused
	for(int i = 0; i < ndims_; i++)
	{
		if (levels_(i)==0)
		{
			indexboundary(i) = 1; // dfixed vector
			nboundary++;
		}
	}
	vector onevec(ndims_); onevec.fill(1.0);
	vector result(ndata); result.fill(0);

	// go through data points
	for(int i = 0; i < ndata; i++)
	{
		vector x = data.row(i);
		vector index = gridutils::SubGridIndex(x, levels_);
		// if no boundary
		if(nboundary == 0)
		{
			int bindex = index.dot(strides_);
			vector relativeX = gridutils::RelativeXScaledSubGrid(x, index, 
																 levels_);
			result[i] += gridutils::HatVal(relativeX)*data_[bindex];
		}
		// else need to visit each boundary basis function
		else
		{
			vector bit(ndims_); bit.fill(0);
			for(int j =0; j < gridutils::PowInt(2,nboundary); j++)
			{
				int bindex = strides_.dot(index+bit);
				vector relativeX = gridutils::RelativeXScaledSubGrid(x, index, 
																	 levels_);
				result[i] += gridutils::HatVal(relativeX)*data_[bindex];
				gridutils::IncreaseBit(bit, onevec, indexboundary);
			}
		}
	}
	return result;

}