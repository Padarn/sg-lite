/*

TEST regular_grid_test: This test will run through basic functionality of the 
		RegularGrid class, making sure that things are initialised correctly, 
		sizes are calculated correctly, and indexing works properly.

*/

#include "sg_lite_include.hpp"
#include "gtest/gtest.h"


// Sets up level and ndim for an anisotropic 
// three dimensional grid.
class RegularGridTest_3D : public ::testing::Test
{
	public:
		int ndim;
		vector level;
		matrix testdata1, testdata2;
		matrix evalpts;

	virtual void SetUp(){
		ndim = 3;
		level.resize(ndim); level.fill(3);
		level(1) = 2;
		testdata1 = read_csv_matrix("data/data3d_one.csv");
		testdata2 = read_csv_matrix("data/data3d_two.csv");
		evalpts.resize(4, ndim);
		evalpts << 0, 0, 0,
				   0, 1, 0,
				   0.6, 0.6, 0.6,
				   0.4, 0.5, 0.4;
	}
};

// Test grid with boundary
TEST_F(RegularGridTest_3D, TestWithBoundary)
{
	RegularGrid grid = RegularGrid(ndim, level, true);

	// ------ Basic test of functionality ------- //

	// test HatDataBin and EvalPoints
	grid.HatDataBin(testdata1);
	vector result = grid.EvalPoints(evalpts);

	EXPECT_EQ(grid.ndata_, 1);
	EXPECT_EQ(result(0), 0);


}

// Test grid without boundary
TEST_F(RegularGridTest_3D, TestWithoutBoundary)
{
	RegularGrid grid = RegularGrid(ndim, level, false);
}

	

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
