#create particle distribution files for testing
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

# Create a cubic grid of particles
# x y z vx vy vz ax ay az rho_0 p_0 
def create_cube_distribution(n_particles_per_side, spacing, filename, rho_0, p_0):
    # Much more efficient: use meshgrid
    coords = np.arange(n_particles_per_side) * spacing
    x, y, z = np.meshgrid(coords, coords, coords, indexing='ij')
    
    # Flatten and stack positions
    positions = np.column_stack([x.ravel(), y.ravel(), z.ravel()])
    
    # Create velocities and accelerations (all zeros)
    n_particles = n_particles_per_side ** 3
    velocities = np.zeros((n_particles, 3))
    accelerations = np.zeros((n_particles, 3))
    
    # Create density and pressure columns
    densities = np.full((n_particles, 1), rho_0)
    pressures = np.full((n_particles, 1), p_0)
    
    # Combine all columns
    particles = np.hstack([positions, velocities, accelerations, densities, pressures])
    
    np.savetxt(filename, particles, delimiter=",", header="x,y,z,vx,vy,vz,ax,ay,az,rho_0,p_0", comments="")

if __name__ == "__main__":
    create_cube_distribution(50, 0.1, r'input/cube_distribution.csv', 100.0, 0.0)
    