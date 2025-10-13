# plot naive example
#read from txt file and plot the data

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.path import Path
import numpy as np
data = np.loadtxt('results/naive_solution.txt')
x = data[:,0]
y = data[:,1]
u = data[:,2]



#no cutoff rectangles
# reshape the area to [0,2] * [-2,2]
x = 2 * (x - np.min(x)) / (np.max(x) - np.min(x))
y = 4 * (y - np.min(y)) / (np.max(y) - np.min(y)) - 2
#no holes
#no Triangulate once
# just plot a rectangle area
tri = mtri.Triangulation(x, y)
# Plot with masked triangulation
plt.figure()
cs = plt.tricontourf(tri, u, levels=128, cmap="RdBu_r")
plt.colorbar(cs)
plt.title("Numerical solution ")
plt.xlabel("x")
plt.ylabel("y")
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig('results/naive_rhs.png', dpi=200, bbox_inches='tight')
plt.show()