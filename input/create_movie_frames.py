import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

#create movie frames for sph simulations
#assumes files named output_particles_step_x.csv where x is the timestep number
# saved in SPH/data/
import os
script_dir = os.path.dirname(os.path.realpath(__file__))
output_dir = os.path.join(script_dir, 'movie_frames')
os.makedirs(output_dir, exist_ok=True)
data_dir = os.path.join(script_dir, '..', 'data')
particle_files = sorted([f for f in os.listdir(data_dir) if f.startswith('output_particles_step_') and f.endswith('.csv')],
                        key=lambda x: int(x.split('_')[-1].split('.')[0]))
for particle_file in particle_files:
    data = np.loadtxt(os.path.join(data_dir, particle_file), delimiter=',', skiprows=1)
    #assume columns are: x,y,z,vx,vy,vz,mass,smoothing_length,density,pressure,sound_speed
    # make a 3d plot with the density ascolor bar

    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]
    density = data[:, 8]
    timestep = particle_file.split('_')[-1].split('.')[0]
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    p = ax.scatter(x, y, z, c=density, cmap='viridis', s=40, vmin = 0, vmax = 1)
    fig.colorbar(p, ax=ax, label='Density')
    ax.set_title(f'SPH Simulation at Timestep {timestep}')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xlim([-5, 5])
    ax.set_ylim([-5, 5])
    ax.set_zlim([-5, 5])
    plt.grid(True)
    frame_filename = os.path.join(output_dir, f'movie_frame_{timestep.zfill(5)}.png')
    plt.savefig(frame_filename, dpi=150)
    plt.close()
print(f'Movie frames saved in directory: {output_dir}')
  
