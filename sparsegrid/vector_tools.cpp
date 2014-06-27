#include "vector_tools.hpp"
#include <iostream>

namespace vectortools{
 
	vector MinVector(vector x, vector y)
	{
	    vector minvec(x.size());
	    for(int i = 0; i < x.size(); i++)
	    {
	    	minvec(i) = fmin(x(i),y(i));
	    }
	    return minvec;
	}

	vector MaxVector(vector x, vector y)
	{
	    vector maxvec(x.size());
	    for(int i = 0; i < x.size(); i++)
	    {
	    	maxvec(i) = fmax(x(i),y(i));
	    }
	    return maxvec;
	}


	int VectorToIndex(vector in, vector maxlevel)
	{
		int index = 0;
		int mult = 1;
		for(int i = 0; i<in.size();i++)
		{
			index += mult*(in(i)+1);
			mult = mult*(maxlevel(i)+2);
		}
		return index;
	}

	vector IndexToVector(int in, vector maxlevel)
	{
		int in2 = in;
		int ndim = maxlevel.size() ;
		vector level(ndim);
		int mult;
		for(int i=0; i<ndim; i++)
		{
			mult = (int) (maxlevel(i)+2);
			level(i) = in % mult - 1;
			in = in / mult;
		}
		return level;
	}

}