#include "sparsegrid.hpp"
#include "utils.cpp"
#include <iostream>


// -------------------------- UTILITIES STRUCT ------------------------------ //

int pow_int(int a,int b){
	return (int) pow(a, b);
}

int boundary_size(const vector & levels, const int &  ndim)
{
	int size = 1;
	for(int i=0; i<ndim; i++){
		size *= pow_int(2,levels(i))+1;
	}
	return size;
}

int no_boundary_size(const vector & levels, const int &  ndim)
{
	int size = 1;
	for(int i=0; i<ndim; i++){
		size *= pow_int(2,levels(i))-1;
	}
	return size;
}

int get_corner_index(vector x, vector levels, bool full, vector & fixedd)
{
	int dim = x.size();
	int index = 0;
	int base;
	int prodd;
	int indexcur;

	// modifier for boundary
	if(full) base = 1;
	else base = -1;

	prodd = 1;

	for(int d = 0; d<dim; d++)
	{
		// If outside domain return -1
		if (x(d)>1 || x(d)<0) return -1;

		indexcur = floor(pow_int(2,levels(d))*x(d));
		// index is one lower if no boundary
		if (!full) indexcur--;

		// deal with no boundary
		// if on left boundary fix and thats all
		if (!full && indexcur ==0){
			fixedd(d)=-1;
		}
		// if one away from right bounday
		// fix and thats all..
		if (!full && indexcur == pow_int(2,levels(d))-2){
			fixedd(d)=1;
			indexcur--;
		}
		// if on right boundary full grid
		if (full && indexcur == pow_int(2,levels(d))){
			fixedd(d)=1;
			indexcur--;
		}
		// right boundary no boundary
		if (!full && indexcur == pow_int(2,levels(d))-1){
			fixedd(d)=-1;
			indexcur--;
		}

		index += prodd*indexcur;
		prodd *= pow_int(2,levels(d))+base;
	}
	return index;

}

vector stride_index(vector levels, bool full)
{
	int base;
	int dprod;
	int d = levels.size();
	vector strides(d);

	// modifier for boundary
	if(full) base = 1;
	else base = -1;

	dprod = 1;
	for(int i = 0; i < d; i++){
		strides(i) = dprod;
		dprod *= pow_int(2,levels(i))+base;
	}

	return strides;

}

void increase_bit(vector & bit, vector zsize){

	int j = 0;
	while(true && j<bit.size())
	{
		if (bit(j)<zsize(j))
		{
			bit(j) += 1;
			break;
		}
		bit(j)=0;
		j++;
	}

}

void increase_bit_exd(vector & bit, vector zsize, int d){

	if (bit.size() == 1)
	{
		return;
	}
	int j = 0;
	while(true)
	{
		if (j == d) j++;
		if (bit(j)<zsize(j))
		{
			bit(j) += 1;
			break;
		}
		bit(j)=0;
		j++;
	}

}

void increase_2_bit_exd(vector & bit, vector fixedd){

	if (bit.size() == 1)
	{
		return;
	}
	int j = 0;
	while(true && j<bit.size())
	{
		if (fixedd(j) != 0) j++;
		else
		{
			if (bit(j)<1)
			{
				bit(j) += 1;
				break;
			}
			bit(j)=0;
			j++;
		}
	}

}

vector get_xleft(vector x, vector levels){

	int n = x.size();
	vector xleft(n);
	for(int i = 0; i < n; i++){
		if (x(i)==1) xleft(i)=1;
		else xleft(i)=pow_int(2,levels(i))*x(i)-floor(pow_int(2,levels(i))*x(i));
	}
	return xleft;

}

vector get_xright_from_xleft(vector xleft){

	int n = xleft.size();
	vector xright(n);
	xright.fill(1);
	xright = xright - xleft;
	return xright;

}

double hatval(vector &bit, vector &xleft, vector &xright)
{
	double val = 1;
	for(int i=0; i<bit.size();i++){
		val*=(1-bit(i))*xright(i) + (bit(i))*xleft(i);
	}
	return val;
}

// NOTE : DOES NOT HANDLE THE SCALE - THIS SHOULD BE COMBINED 
// AT THE END;
void tridiagonal_projection(vector & beta, int index_start, int index_stride,
						    int size){
	

	// size
	int N = size;
	
	// Modify inplace beta to do projection (has basis) for tridiagonal
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

	for(int i = N-2; i >= 0; i--)
	{
		beta(zero+i*step) -= c(i)*beta(zero+(i+1)*step);
	}

}
// -------------------------------------------------------------------------- //


// ------------------------- REGULARGRID STRUCT ----------------------------- //

// Constructor
RegularGrid::RegularGrid(int ndims, vector levels, bool boundary)
{

	ndims_ = ndims;
	levels_ = levels;
	boundary_ = boundary;


}

// Initialize the data array
void RegularGrid::Initialize()
{
	// NOTE: level 0 means boundary
	//       level 1 add central point
	//       level d refine
	// Thus no boundary and level 0 is degenerate
	if(boundary_)
	{
		size_ = boundary_size(levels_,ndims_);
	}
	else
	{
		size_ = no_boundary_size(levels_,ndims_);
	}

	data_ = vector(size_);
	data_.fill(0);
}

// Evaluate data on regular grid hat functions
void RegularGrid::EvaluateData(const matrix & data)
{
	
	// Set data to zero (intialize?)
	data_.fill(0);
	int Ndata = 0;

	// Set up data loop
	int ndata = data.rows();
	vector strides = stride_index(levels_, boundary_);

	for(int i = 0; i < ndata; i++)
	{
		vector x = data.row(i);
		vector fixedd(ndims_);
		fixedd.fill(0);
		int indexbase = get_corner_index(data.row(i), levels_, boundary_,
										 fixedd);
		// deal with out sample data
		if (indexbase > -1 and indexbase < size_)
		{
			Ndata++;
			vector xleft = get_xleft(x, levels_);
			vector xright =  get_xright_from_xleft(xleft);
			//hat_eval_right(xright);
			//hat_eval_left(xleft);


			// go through all dimensions unrolled - avoid recursion
			vector bit(ndims_);
			bit.fill(0);
			// set fixedd in place - fixedd stores 1 for a fixed right and -1
			// for a fixed left
			int nfixed = 0;
			for(int j = 0; j< ndims_;j++)
			{
				if (fixedd(j) == -1){
					bit(j)=0;
					nfixed++;
				}else if(fixedd(j)==1){
					nfixed++;
					bit(j)=1;
				}
			}

			for(int j = 0; j<pow_int(2,ndims_- nfixed); j++)
			{
				int index; 
				index = indexbase + bit.dot(strides);
				double val = hatval(bit,xleft,xright);
				increase_2_bit_exd(bit,fixedd);	
				data_[index] += val;

			}
		}
	}
	data_ = data_*1.0/Ndata;
}

// Do L2 projection of stored data
void RegularGrid::ProjectData()
{	
	// Get vector of strides
	vector strides = stride_index(levels_, boundary_);
	vector zsize(ndims_);
	for( int d=0; d<ndims_; d++){
		if (boundary_) zsize(d) = pow_int(2,levels_(d))+1;
		else zsize(d) = pow_int(2,levels_(d))-1;
	}

	// For each dimension
	for(int d=0; d<ndims_; d++)
	{
		// Interate through 'start' indicies
		int dsize = zsize(d);
		int remsize = size_/dsize;

		int stride = strides(d);
		vector bit(ndims_);
		bit.fill(0);

		for(int i=0; i < remsize; i++)
		{
			//get indexstart from bit method
			int indexstart = bit.dot(strides);
			increase_bit_exd(bit, zsize, d);
			// Do projection based at index of length specified.
			tridiagonal_projection(data_, indexstart, stride, dsize);
		}

		// Scale
		data_ = data_*pow_int(2, levels_(d)-1);
	}

}

// Calculate value of grid at any point
vector RegularGrid::EvalPoints(matrix data)
{
	
	vector result(data.rows());
	result.fill(0);
	// Set up data loop
	int ndata = data.rows();
	vector strides = stride_index(levels_, boundary_);

	for(int i = 0; i < ndata; i++)
	{
		vector x = data.row(i);
		vector fixedd(ndims_);
		fixedd.fill(0);
		int indexbase = get_corner_index(data.row(i), levels_, boundary_,
										 fixedd);

		if (indexbase > -1 and indexbase < size_)
		{
			vector xleft = get_xleft(x, levels_);
			vector xright =  get_xright_from_xleft(xleft);
			//hat_eval_right(xright);
			//hat_eval_left(xleft);


			// go through all dimensions unrolled - avoid recursion
			vector bit(ndims_);
			bit.fill(0);
			// set fixedd in place - fixedd stores 1 for a fixed right and -1
			// for a fixed left
			int nfixed = 0;
			for(int j = 0; j< ndims_;j++)
			{
				if (fixedd(j) == -1){
					bit(j)=0;
					nfixed++;
				}else if(fixedd(j)==1){
					nfixed++;
					bit(j)=1;
				}
			}

			for(int j = 0; j<pow_int(2,ndims_- nfixed); j++)
			{
				int index; 
				index = indexbase + bit.dot(strides);
				double val = hatval(bit,xleft,xright);
				result(i) += val*data_(index);
				increase_2_bit_exd(bit,fixedd);
			}
		} 

	}
	return result;
}

// Calculate value of grid at any point
vector RegularGrid::EvalPointsGrid(int res)
{
	// make grid of data
	int N = pow(res,ndims_);
	matrix grid(N, ndims_);
	vector dx(ndims_);
	dx.fill(1.0/(res-1));
	vector bit(ndims_);
	bit.fill(0);
	vector len(ndims_);
	len.fill(res-1);
	grid.fill(0);
	vector x(ndims_);

	for(int i = 0; i < N; i++){
		vector x = bit.cwiseProduct(dx);
		grid.row(i)=x;
		increase_bit(bit, len);
	}

	return EvalPoints(grid);
}

// -------------------------------------------------------------------------- //



// --------------------------- SUBGRID STRUCT ------------------------------- //





// -------------------------------------------------------------------------- //


// ----------------------- COMBINATIONGRID STRUCT --------------------------- //

int fact(int n)
{
	int fn = 1;
	for(int i = 1;i<=n;i++)
	{
		fn = fn * i;
	}
	return fn;
}

int nchoosek(int n, int k){
	int p1 = fact(n);
	int p2 = fact(n-k);
	int p3 = fact(k);
	int result = p1/(p2*p3);
	return result;
}


int findnonzero(vector vecin, int num)
{
	int n = vecin.size();
	int found = 0;
	int i = n-1;
	while(i>=0)
	{
		if (vecin(i)!=0)
		{
			found++;
			if (found == num)
			{
				return i;
			}
		}
		i--;
	}
	return -1;
}

void advance_levels(vector & levels)
{
	int n = levels.size();
	if (n==1) return;

	int lastnonzero = findnonzero(levels,1);
	if (lastnonzero < n-1)
	{
		levels(lastnonzero) = levels(lastnonzero)-1;
		levels(lastnonzero+1) = 1;
	} else 
	{
		int indnonzero = findnonzero(levels,2);
		int delta = levels(lastnonzero);
		levels(lastnonzero)=0;
		int newval = delta + 1;
		levels(indnonzero) = levels(indnonzero)-1;
		levels(indnonzero+1) = newval; 
	}
}

CombinationGrid::CombinationGrid(int ndims, bool boundary)
{
	ndims_ = ndims;
	boundary_ = boundary;
}


// ---------------------
//     BUILD GRIDS
// ---------------------

void CombinationGrid::SetupStandardCombinationGrid(int maxlevel)
{
	// NOTE: Does not intialize the grids. Only do this if they
	// are needed.

	/// NOTE :: CURRENT ASSUMING BOUNDARY

	// TO CHANGE TO NON - BOUNDARY - NEED TO CHANGE LEVELS TO 
	// START COUNTING AT 1... or just minus 1 to top level and
	// add one to level?

	int coef;
	vector level(ndims_);
	vector ones(ndims_);
	if (!boundary_) ones.fill(1);
	else ones.fill(0);

	if(!boundary_) maxlevel--;

	for(int curlvl = maxlevel - ndims_ + 1; curlvl <= maxlevel; curlvl++)
	{
		coef = pow(-1,ndims_-curlvl)*nchoosek(ndims_-1,maxlevel-curlvl);
		level.fill(0);
		level(0)=curlvl;
		int n_atlevel = nchoosek(ndims_+curlvl-1,curlvl);
		//vector newlevel = level+ones; // CHECK IF NEEDED
		levels_.push_back(level+ones);
		coefs_.push_back(coef);
		for(int j=0; j<n_atlevel-1; j++){
			advance_levels(level);
			//newlevel = level+ones; // CHECK IF NEEDED
			levels_.push_back(level+ones);
			coefs_.push_back(coef);
		}
	}

}

void CombinationGrid::Initialize()
{

	// Go through all levels and make a 'regular grid'
	for(unsigned int i = 0; i < levels_.size(); i++)
	{
		regulargridptr grid(new RegularGrid(ndims_, levels_[i], boundary_));
		grid->Initialize();
		grids_.push_back(grid);
	}

}

// ---------------------
//   DO STUFF TO GRIDS
// ---------------------

// Evaluate data on grid
void CombinationGrid::EvaluateData(const matrix & data)
{
	for(unsigned int i = 0; i< grids_.size(); i++)
	{
		grids_[i]->EvaluateData(data);
	}
}

// L2 Projection of data (hat functions)
void CombinationGrid::ProjectData()
{
	for(unsigned int i = 0; i< grids_.size(); i++)
	{
		grids_[i]->ProjectData();
	}
}

// ---------------------
//     EVAL GRIDS
// ---------------------

vector CombinationGrid::EvalPoints(matrix data)
{
	vector result(data.rows());
	result.fill(0);
	for(unsigned int i = 0; i< grids_.size(); i++)
	{
		vector part = grids_[i]->EvalPoints(data);
		result = result + coefs_[i]*part;
	}
	return result;
}

vector CombinationGrid::EvalPointsGrid(int res)
{
	// TODO MAKE THIS MORE ROBUST
	vector part = grids_[0]->EvalPointsGrid(res);
	vector result = coefs_[0]*part;
	for(unsigned int i = 1; i< grids_.size(); i++)
	{
		vector part = grids_[i]->EvalPointsGrid(res);
		result = result + coefs_[i]*part;
	}
	return result;
}



// -------------------------------------------------------------------------- //

