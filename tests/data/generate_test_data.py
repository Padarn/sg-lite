"""
    generate_test_data.py
    
    Generates testing data from a specified seed. This allows the same test
    data to be used without storing the data.

    Padarn Wilson : 25/06/14
"""


import numpy as np
np.random.seed(1234)

# make first 3d test data - just one points at
# (0.6,0.6,0.6)
data3d = np.array([[0.6, 0.6, 0.6]])
np.savetxt('data3d_one.csv', data3d, delimiter=',')

# make second 3d test data - three points in corners
# (0,0,0), (1,0,0), (0,0.5,0), (1,1,1)
data3d = np.array([[0.0, 0.0, 0.0],
                   [1.0, 0.0, 0.0],
                   [0.0, 0.5, 0.0],
                   [1.0, 1.0, 1.0]])
np.savetxt('data3d_two.csv', data3d, delimiter=',')
