#read from txt file and plot the data
# use matplotlib to plot the data in results/solution.txt
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.path import Path
import numpy as np
data = np.loadtxt('results/solution.txt')
x = data[:,0]
y = data[:,1]
u = data[:,2]
# reshape the area to [0,2] * [-2,2]
x = 2 * (x - np.min(x)) / (np.max(x) - np.min(x))
y = 4 * (y - np.min(y)) / (np.max(y) - np.min(y)) - 2

# Define hole polygons (your three removed triangles)
holes = [
    np.array([[0.0, 1.0], [0.0, 2.0], [1.0, 2.0]]),   # T1
    np.array([[1.0, 2.0], [2.0, 2.0], [2.0, 1.0]]),   # T2
    np.array([[1.0,-1.0], [2.0, 1.0], [2.0,-2.0]]),   # T3
]
hole_paths = [Path(h) for h in holes]

# Triangulate once
tri = mtri.Triangulation(x, y)

# Mask triangles whose centroid falls inside any hole
tri_pts_x = x[tri.triangles]
tri_pts_y = y[tri.triangles]
xc = tri_pts_x.mean(axis=1)
yc = tri_pts_y.mean(axis=1)

mask = np.zeros(tri.triangles.shape[0], dtype=bool)
for hp in hole_paths:
    mask |= hp.contains_points(np.c_[xc, yc])

tri.set_mask(mask)

# Plot with masked triangulation
plt.figure()
cs = plt.tricontourf(tri, u, levels=128, cmap="RdBu_r")
plt.colorbar(cs)
plt.title("Numerical solution")
plt.xlabel("x")
plt.ylabel("y")
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig('results/solution.png', dpi=200, bbox_inches='tight')
plt.show()