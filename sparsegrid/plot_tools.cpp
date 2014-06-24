#include "plot_tools.hpp"
#include <iostream>

namespace plottools{

	matrix MeshGrid(int res, int d)
	{

		int N = gridutils::PowInt(res, d);
		matrix grid(N, d); grid.fill(0);
		vector bit(d); bit.fill(0);
		vector len(d); len.fill(res-1);
		vector dx(d); dx.fill(1.0/(res-1));
		vector x(d);

		for(int i = 0; i < N; i++){
			vector x = bit.cwiseProduct(dx);
			grid.row(i)=x;
			gridutils::IncreaseBit(bit, len, -1);
		}

		return grid;
	}
	
}