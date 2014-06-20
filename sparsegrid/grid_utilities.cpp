#include "grid_utilities.hpp"
#include <iostream>

namespace gridutils {

	// ----------------------------------------------            
	//               GENERIC UTILITIES
	// ----------------------------------------------

	int PowInt(int a, int b)
	{
		int result = 1;
		for(int i = 0; i<b; i++){
			result = result*b;
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

	bool BoundaryBitAndFixed(vector strides, vector corner, bool boundary, 
							 vector & bit, vector & fixed, int & nfixed)
	{
		int d = strides.size();
		for(int i = 0; i<d; i++)
		{
			int dindex = corner(i);
			if (dindex < 0)
			{
				if (dindex == -1 && !boundary)
				{
					bit(i) = 1;
					fixed(i) = 1;
					nfixed++;
				}
				else return false;
			}
			if (dindex>=strides(i))
			{
				if(dindex == strides(i))
				{
					bit(i) = 0;
					fixed(i) = 1;
					nfixed++;
				}
			else return false;
			} 
		}

		return true;
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
	// ----------------------------------------------            
	//                  BASIS EVAL
	// ----------------------------------------------

	double HatVal(vector x)
	{
		double val = 1;
		int d = x.size();
		for(int i=0; i<d;i++){
			val*=(1-abs(x(i)));
		}
		return val;
	}


}