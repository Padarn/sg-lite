import numpy as np
import matplotlib.pyplot as plt
import sys
from mpl_toolkits.mplot3d import Axes3D


filename = sys.argv[1]
data = np.loadtxt(filename)

n = np.int(np.sqrt(len(data)))
data = data.reshape(n,n)

X, Y = np.meshgrid(range(n),range(n))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(X,Y,data, rstride=1, cstride=1, alpha = 0.25)

plt.show()