import numpy as np
import matplotlib.pyplot as plt
import h5py

# Load h5 file
dist = h5py.File('cube_distribution.h5', 'r')
positions = dist['positions'][:]
masses = dist['mass'][:]
spacing = (positions.max(axis=0) - positions.min(axis=0)).max() / (len(np.unique(positions[:,0])) - 1)
density = masses / spacing**3  # assuming uniform distribution for simplicity 
# Plot distribution with density as heatmap
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
sc = ax.scatter(positions[:,0], positions[:,1], positions[:,2], c=density, s=40, cmap='viridis')
ax.set_title(f'Particle Distribution')
plt.colorbar(sc, label='Density')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.savefig('particle_distribution.png', dpi=300)
plt.show()