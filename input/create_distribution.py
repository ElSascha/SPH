#create particle distribution files for testing
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import h5py

# Create a cubic grid of particles saved in a h5 file that rotates around the z axis
# x y z vx vy vz mass smoothing_length  
def create_cube_distribution(n_particles_per_side, size, filename, rho_0):
    # Much more efficient: use meshgrid
    spacing = size / n_particles_per_side
    coords = np.arange(n_particles_per_side) * spacing
    x, y, z = np.meshgrid(coords, coords, coords, indexing='ij')
    
    # Flatten and stack positions
    positions = np.column_stack([x.ravel(), y.ravel(), z.ravel()])
    
    # Create velocities for rotation around z axis
    n_particles = n_particles_per_side ** 3
    velocities = np.zeros((n_particles, 3))
    
    # Center should be the actual center of the particle distribution, not size/2
    center_x = coords.mean()
    center_y = coords.mean()
    
    omega = 1.0  # Drehgeschwindigkeit
    for i in range(n_particles):
        dx = positions[i,0] - center_x
        dy = positions[i,1] - center_y
        velocities[i,0] = -omega * dy  # vx = -ω y
        velocities[i,1] = omega * dx   # vy = ω x
        velocities[i,2] = 0.0
    
    # Verify zero net momentum
    print(f"Net momentum: vx={velocities[:,0].mean():.2e}, vy={velocities[:,1].mean():.2e}")
    
    # Create density and pressure columns
    mass = (rho_0 * spacing**3)
    smoothing_length = 1.2 * spacing  # typical choice
    
    with h5py.File(filename, 'w') as f:
        f.create_dataset('positions', data=positions)
        f.create_dataset('velocities', data=velocities)
        f.create_dataset('mass', data=mass * np.ones(n_particles))
        f.create_dataset('smoothing_length', data=smoothing_length * np.ones(n_particles))
    #np.savetxt(filename, particles, delimiter=",", header="x,y,z,vx,vy,vz,mass,smoothing_length", comments="")
    print(f"Created cube distribution with {n_particles} particles in '{filename}'")

import numpy as np
import h5py

def shear_patch_test_2d(n, size, shear_rate, rho_0 ,filename):
    """
    Creates a 2D shear patch test configured exactly like in SPH solid mechanics papers.
    
    n           number of particles per side
    size        physical size of the square domain
    shear_rate  gamma (in v_x = gamma * y)
    filename    output .h5
    """
    # 2D grid
    spacing = size / (n - 1)
    coords = np.arange(n) * spacing
    x, y = np.meshgrid(coords, coords, indexing='ij')
    
    positions = np.column_stack([x.ravel(), y.ravel(), np.zeros(n*n)])
    
    velocities = np.zeros_like(positions)
    velocities[:,0] = shear_rate * positions[:,1]   # vx = gamma * y
    
    # Typical SPH choices
    mass = rho_0 * spacing**2
    smoothing_length = 1.2 * spacing
    
    with h5py.File(filename, 'w') as f:
        f.create_dataset('positions', data=positions)
        f.create_dataset('velocities', data=velocities)
        f.create_dataset('mass', data=mass * np.ones(n*n))
        f.create_dataset('smoothing_length', data=smoothing_length * np.ones(n*n))
    
    print(f"[OK] 2D shear patch test generated with {n*n} particles → {filename}")

if __name__ == "__main__":
    create_cube_distribution(6, 5.0, r'cube_distribution.h5', 2900.0)
    #shear_patch_test_2d(30, 1.0, 1.0, 1.0, r'shear_patch_test.h5')
