import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# Example data: Replace this with your list of coordinates
coordinates = [
    [1, 2, 3],
    [1, 3, 4],
    [2, 2, 5],
    [2, 3, 6],
    [3, 2, 7],
    [3, 3, 8],
    [3, 4, 9],
]

# Convert the list to a NumPy array
points = np.array(coordinates)

# Separate x, y, z values
x = points[:, 0]
y = points[:, 1]
z = points[:, 2]

# Create a grid for x and y
x_grid = np.linspace(min(x), max(x), 50)
y_grid = np.linspace(min(y), max(y), 50)
X, Y = np.meshgrid(x_grid, y_grid)

# Interpolate z values over the grid
Z = griddata((x, y), z, (X, Y), method='linear')

# Plot the surface
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='k')

# Add labels and color bar
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Z-axis')
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.title('3D Surface Plot from Coordinates')
plt.show()
