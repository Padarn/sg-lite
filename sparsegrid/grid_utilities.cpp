#include "grid_utilities.hpp"
#include <iostream>
#include <cmath>

namespace gridutils {

	// ----------------------------------------------            
	//               GENERIC UTILITIES
	// ----------------------------------------------

	int PowInt(int a, int b)
	{
		int result = 1;
		for(int i = 0; i<b; i++){
			result = result*a;
		}
		return result;
	}

	// ----------------------------------------------            
	//                  GRID SIZES 
	// ----------------------------------------------

	int RegularGridSize1D(int level, bool boundary)
	{
		if (boundary){
			return PowInt(2,level)+1;
		} else {
			return PowInt(2,level)-1;
		}
	}

	int SubGridSize1D(int level)
	{
		if (level == 0) return 2;
		else return PowInt(2,level-1);
	}

	int RegularGridSize(vector level, bool boundary)
	{
		int n = level.size();
		int size = 1;
		for(int i = 0; i<n; i++)
		{
			size *= RegularGridSize1D(level(i), boundary);
		}
		return size;
	}

	int SubGridSize(vector level)
	{
		int n = level.size();
		int size = 1;
		for(int i = 0; i<n; i++)
		{
			size *= SubGridSize1D(level(i));
		}
		return size;
	}

	vector RegularGridSizeVector(vector level, bool boundary)
	{
		int n = level.size();
		vector sizes(n);
		for(int i = 0; i<n; i++)
		{
			sizes(i)=RegularGridSize1D(level(i), boundary);
		}
		return sizes;
	}

	vector SubGridSizeVector(vector level)
	{
		int n = level.size();
		vector sizes(n);
		for(int i = 0; i<n; i++)
		{
			sizes(i)=SubGridSize1D(level(i));
		}
		return sizes;
	}

	vector StrideFromSize(vector sizes)
	{
		vector strides(sizes.size());
		strides(0) = 1;
		for(int i = 0; i<sizes.size()-1; i++)
		{
			strides(i+1)=sizes(i);
		}
		return strides;
	}

	// ----------------------------------------------            
	//               INDEX GENERATION
	// ----------------------------------------------

	void IncreaseBit(vector & bit, vector dimsizes, int d)
	{
		int j = 0;
		while(j<bit.size())
		{
			if (j == d) j++;
			else
			{
				if (bit(j)<dimsizes(j))
				{
					bit(j) += 1;
					break;
				}
				bit(j)=0;
				j++;
			}
		}
	}

	void IncreaseBit(vector & bit, vector dimsizes, vector dvec)
	{
		int j = 0;
		int i = 0; // keeps track of held dimensions
		while(j<bit.size())
		{
			if (dvec(i)){i++; j++;} // d 0 for free, 1 for fixed
			else
			{
				if (bit(j)<dimsizes(j))
				{
					bit(j) += 1;
					break;
				}
				bit(j)=0;
				j++;
				i++;
			}
		}
	}

	void DistributeBit(vector & bit, vector sizes)
	{
		if (bit.sum() > sizes.sum())
		{
			bit.fill(0);
			return;
		}
		for(int i = 0; i < sizes.size()-1; i++)
		{
			if(bit(i) > sizes(i))
			{
				int diff = bit(i)-sizes(i);
				bit(i+1) = bit(i+1) + diff;
				bit(i) = sizes(i);
			}
		}
	}

	void IncreaseBitMonotone(vector & bit, vector sizes)
	{
		int n = bit.size();
		while(true)
		{
			// Get non-zero entires
			int last = -1;
			int second = -1;
			int numnonzero=0;
			for(int i = 0; i < n; i++)
			{
				if(bit(i)!=0)
				{
					second = last;
					last = i;
					numnonzero++;
				}
			}

			// If the zero vector
			if (numnonzero==0){
				bit(0)=1;
				return;
			}

			// If non-zero is not on the end move over 
			else if (last != n-1)
			{
				bit(last) = bit(last) - 1;
				bit(last + 1) = 1;
				DistributeBit(bit, sizes);
			}

			// If is on end and only one
			else if (last == n-1 && numnonzero == 1)
			{
				bit(0) = bit(n-1) + 1;
				bit(n-1) = 0;
				DistributeBit(bit, sizes);
			}

			// Otherwise
			else
			{
				int tmp = bit(last);
				bit(last) = 0;
				bit(second) = bit(second) - 1;
				bit(second + 1) = bit(second + 1) + tmp + 1;
				DistributeBit(bit, sizes);
			}

			if (!(bit(n-1) > sizes(n-1)))
			{
				return;
			}
		}
	}

	vector SubGridIndex(vector x, vector level)
	{
		int d = level.size();
		vector index(d);
		for(int i=0; i<d; i++)
		{
			if (level(i)==0) index(i)=0;
			else index(i)=floor(PowInt(2,level(i)-1) * x(i));
		}
		return index;
	}

	vector CornerStrides(vector x, vector level, bool boundary)
	{
		int d = level.size();
		// Initialize result
		vector corner(d);
		for(int i=0; i<d; i++)
		{
			corner(i) = floor(PowInt(2, level(i)) * x(i));
			// (1-boundary) shift accounts for boundary index change
			corner(i) = corner(i) - (1-boundary);
		}
		return corner;
	}

	// TODO - NEEDS FIXING
	bool BoundaryBitAndFixed(vector strides, vector corner, bool boundary, 
							 vector & bit, vector & fixed, int & nfixed)
	{
		bool indomain = true;
		int d = strides.size();
		for(int i = 0; i<d; i++)
		{
			int dindex = corner(i);
			if (dindex < 0)
			{
				if (dindex == -1)
				{
					bit(i) = 1;
					fixed(i) = 1;
					nfixed++;
				}
				if (boundary) indomain = false;
			}
			if (dindex>=(strides(i)-1))
			{
				if(dindex == (strides(i)-1))
				{
					bit(i) = 0;
					fixed(i) = 1;
					nfixed++;
				}
			else indomain = false;
			} 
		}

		return indomain;
	}

	// ----------------------------------------------            
	//                  X POSITION
	// ----------------------------------------------

	vector IndexX(vector index, vector level, bool boundary)
	{
		int d = level.size();
		// Initialize result
		vector X(d);
		for(int i=0; i<d; i++)
		{
			X(i) = index(i) + (1-boundary);
			X(i) = X(i)*pow(2,-level(i));
		}
		return X;

	}

	vector IndexXSubGrid(vector index, vector level)
	{
		int d = level.size();
		// Initialize result
		vector X(d);
		for(int i=0; i<d; i++)
		{
			if (level(i)==0){
				X(i)=index(i);
			} else {
				X(i)=(2*index(i)+1)*pow(2,-level(i));
			}
		}
		return X;

	}

	vector RelativeX(vector x, vector index, vector level, bool boundary)
	{
		int d = level.size();
		// Initialize result
		vector relativeX(d);
		vector indexX = IndexX(index, level, boundary);
		for(int i=0; i<d; i++)
		{
			relativeX(i) = x(i) - indexX(i);
		}
		return relativeX;
	}

	vector RelativeXScaled(vector x, vector index, vector level, bool boundary)
	{	
		int d = level.size();
		// Initialize result
		vector relativeX(d);
		vector indexX = IndexX(index, level, boundary);
		for(int i=0; i<d; i++)
		{
			relativeX(i) = (x(i) - indexX(i))*PowInt(2, level(i));
		}
		return relativeX;	
	}

	vector RelativeXScaledSubGrid(vector x, vector index, vector level)
	{	
		int d = level.size();
		// Initialize result
		vector relativeX(d);
		vector indexX = IndexXSubGrid(index, level);
		for(int i=0; i<d; i++)
		{
			relativeX(i) = (x(i) - indexX(i))*PowInt(2, level(i));	
		}
		return relativeX;	
	}

	// ----------------------------------------------            
	//                  BASIS EVAL
	// ----------------------------------------------

	double HatVal(vector x)
	{
		double val = 1;
		int d = x.size();
		for(int i=0; i<d;i++){
			val*=(1-fabs(x(i)));
		}
		return val;
	}

	double HaarVal(vector x)
	{
		double val = 1;
		int d = x.size();
		for(int i=0; i<d;i++){
			if (x(i)>0) val*=(-1);
		}
		return val;	
	}


}