#create particle distribution files for testing
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import h5py

# Create a cubic grid of particles saved in a h5 file
# x y z vx vy vz mass smoothing_length  
def create_cube_distribution(n_particles_per_side, size, filename, rho_0):
    # Much more efficient: use meshgrid
    spacing = size / n_particles_per_side
    coords = np.arange(n_particles_per_side) * spacing
    x, y, z = np.meshgrid(coords, coords, coords, indexing='ij')
    
    # Flatten and stack positions
    positions = np.column_stack([x.ravel(), y.ravel(), z.ravel()])
    
    # Create velocities and accelerations (all zeros)
    n_particles = n_particles_per_side ** 3
    velocities = np.zeros((n_particles, 3))
    
    # Create density and pressure columns
    mass = (rho_0 * spacing**3)
    smoothing_length = 1.2 *spacing  # typical choice
    
    with h5py.File(filename, 'w') as f:
        f.create_dataset('positions', data=positions)
        f.create_dataset('velocities', data=velocities)
        f.create_dataset('mass', data=mass * np.ones(n_particles))
        f.create_dataset('smoothing_length', data=smoothing_length * np.ones(n_particles))
    #np.savetxt(filename, particles, delimiter=",", header="x,y,z,vx,vy,vz,mass,smoothing_length", comments="")
    print(f"Created cube distribution with {n_particles} particles in '{filename}'")

if __name__ == "__main__":
    create_cube_distribution(10, 5.0, r'cube_distribution.h5', 2900.0)

