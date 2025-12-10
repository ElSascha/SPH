import numpy as np
import os
import matplotlib.pyplot as plt

script_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(script_dir, '..', 'data')
particle_files = sorted([f for f in os.listdir(data_dir) if f.startswith('output_particles_step_') and f.endswith('.csv')],
                        key=lambda x: int(x.split('_')[-1].split('.')[0]))

timesteps = []
angular_momenta = []

for particle_file in particle_files:
    data = np.loadtxt(os.path.join(data_dir, particle_file), delimiter=',', skiprows=1)
    x, y, z = data[:, 0], data[:, 1], data[:, 2]
    vx, vy, vz = data[:, 3], data[:, 4], data[:, 5]
    mass = data[:, 6]
    
    # Compute center of mass
    total_mass = np.sum(mass)
    com_x = np.sum(mass * x) / total_mass
    com_y = np.sum(mass * y) / total_mass
    com_z = np.sum(mass * z) / total_mass
    
    # Positions relative to center of mass
    x_rel = x - com_x
    y_rel = y - com_y
    z_rel = z - com_z
    
    # Compute angular momentum about center of mass: L = sum(m * r_rel x v)
    Lx = np.sum(mass * (y_rel * vz - z_rel * vy))
    Ly = np.sum(mass * (z_rel * vx - x_rel * vz))
    Lz = np.sum(mass * (x_rel * vy - y_rel * vx))
    
    timestep = int(particle_file.split('_')[-1].split('.')[0])
    timesteps.append(timestep)
    angular_momenta.append(np.array([Lx, Ly, Lz]))

initial_L = angular_momenta[0]
initial_L_mag = np.linalg.norm(initial_L)
print(f'Initial Angular Momentum (about CoM): Lx={initial_L[0]:.6e}, Ly={initial_L[1]:.6e}, Lz={initial_L[2]:.6e}, |L|={initial_L_mag:.6e}')

relative_changes = []
for i, (t, L) in enumerate(zip(timesteps, angular_momenta)):
    delta_L = np.linalg.norm(L - initial_L)
    rel_change = delta_L / initial_L_mag if initial_L_mag != 0 else 0
    relative_changes.append(rel_change)
    if i > 0:
        print(f'Timestep {t}: ΔL magnitude={delta_L:.6e}, Relative Change={rel_change:.6e}')

print('Angular momentum check complete.')

# Plot
plt.figure(figsize=(8, 5))
plt.plot(timesteps, relative_changes, marker='o')
plt.xlabel('Timestep')
plt.ylabel('Relative Change in |L|')
plt.title('Angular Momentum Conservation Over Time (about CoM)')
plt.grid()
plt.savefig('angular_momentum_conservation.png', dpi=300)
plt.show()