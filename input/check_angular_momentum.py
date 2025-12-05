import numpy as np
#check angular momentum conservation in sph simulation outputs for every timestep
import os
script_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(script_dir, '..', 'data')
particle_files = sorted([f for f in os.listdir(data_dir) if f.startswith('output_particles_step_') and f.endswith('.csv')],
                        key=lambda x: int(x.split('_')[-1].split('.')[0]))
initial_angular_momentum = None
for particle_file in particle_files:
    data = np.loadtxt(os.path.join(data_dir, particle_file), delimiter=',', skiprows=1)
    #assume columns are: x,y,z,vx,vy,vz,mass,smoothing_length,density,pressure,sound_speed
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]
    vx = data[:, 3]
    vy = data[:, 4]
    vz = data[:, 5]
    mass = data[:, 6]
    #compute angular momentum L = r x p = r x (m*v)
    Lx = np.sum(mass * (y * vz - z * vy))
    Ly = np.sum(mass * (z * vx - x * vz))
    Lz = np.sum(mass * (x * vy - y * vx))
    L = np.array([Lx, Ly, Lz])
    if initial_angular_momentum is None:
        initial_angular_momentum = L
        print(f'Initial Angular Momentum: Lx={Lx}, Ly={Ly}, Lz={Lz}')
    else:
        delta_L = L - initial_angular_momentum
        delta_L_magnitude = np.linalg.norm(delta_L)
        initial_L_magnitude = np.linalg.norm(initial_angular_momentum)
        relative_change = delta_L_magnitude / initial_L_magnitude if initial_L_magnitude != 0 else 0
        timestep = particle_file.split('_')[-1].split('.')[0]
        print(f'Timestep {timestep}: ΔL magnitude={delta_L_magnitude:.6e}, Relative Change={relative_change:.6e}')
print('Angular momentum check complete.')