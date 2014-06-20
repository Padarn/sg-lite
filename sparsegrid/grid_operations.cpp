#include "grid_operations.hpp"

namespace gridops{

	void Projection1DSlice(vector & beta, int index_start, int index_stride,
						    int size, int level)
	{	
		int N = size;	
		double cp0 = 1.0/2.0;
		vector c(N-1);
		c.fill(1.0/6);
		c(0)=cp0;

		int zero = index_start;
		int step = index_stride;
		beta(zero) = beta(zero)*3;

		for(int i = 1; i < N -1; i++)
		{
			c(i) = c(i)/(2.0/3-1.0/6*c(i-1));
			beta(zero+i*step)=(beta(zero+i*step)-1.0/6.0*beta(zero+(i-1)*step))/
							  (2.0/3.0-1.0/6.0*c(i-1));
		}

		beta(zero+(N-1)*step)=(beta(zero+(N-1)*step)-1.0/6.0*beta(zero+((N-1)-1)*step))/
							  (1.0/3.0-1.0/6.0*c((N-1)-1));

		// multiply by scale and it will propagate down
		beta(zero+(N-1)*step) *= gridutils::PowInt(2,level-1);
		for(int i = N-2; i >= 0; i--)
		{
			beta(zero+i*step) -= c(i)*beta(zero+(i+1)*step);
		}
	}

	void Projection1D(vector & data, vector strides, vector levels, int size,
					int dim)
	{
		int ndim = levels.size();
		int dsize;
		int stride;

		// Interate through 'start' indicies and project
		vector zsize = strides;
		zsize.array()-1;
		dsize = size*1.0/strides(dim);
		vector bit(ndim); bit.fill(0);
		stride = strides(dim);

		for(int i=0; i < dsize; i++)
		{
			int indexstart = bit.dot(strides);
			Projection1DSlice(data, indexstart, stride, dsize, levels(dim));
			gridutils::IncreaseBit(bit, zsize, dim);
		}

	}

	void ProjectionND(vector & data, vector strides, vector levels, int size)
	{

		int ndim = levels.size();
		for(int d=0; d<ndim; d++)
		{
			Projection1D(data, strides, levels, size, d);
		}
	}
}