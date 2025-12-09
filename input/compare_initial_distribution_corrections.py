import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data_raw = pd.read_csv('../output_particles_initial.csv')
data = pd.read_csv('../output_particles_S.csv')
data_cs = pd.read_csv('../output_particles_CS.csv')

# Set the color bar range from 50 to 100
vmin = 1800
vmax = 2940

# make 3 3d scatter plots of the initial distribution next to each other with density as color map
fig = plt.figure(figsize=(18, 6))   
ax1 = fig.add_subplot(131, projection='3d')
sc1 = ax1.scatter(data_raw['x'], data_raw['y'], data_raw['z'], c=data_raw['density'], s=40, cmap='viridis', vmin=vmin, vmax=vmax)
ax1.set_title('Initial Distribution (Uncorrected Data)')
plt.colorbar(sc1, ax=ax1, label='Density')
ax1.set_xlabel('X')
ax1.set_ylabel('Y')
ax1.set_zlabel('Z')

ax2 = fig.add_subplot(132, projection='3d')
sc2 = ax2.scatter(data['x'], data['y'], data['z'], c=data['density'], s=40, cmap='viridis', vmin=vmin, vmax=vmax)
ax2.set_title('Initial Distribution (Shepard Corrected Data)')
plt.colorbar(sc2, ax=ax2, label='Density')
ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.set_zlabel('Z')

ax3 = fig.add_subplot(133, projection='3d')
sc3 = ax3.scatter(data_cs['x'], data_cs['y'], data_cs['z'], c=data_cs['density'], s=40, cmap='viridis', vmin=vmin, vmax=vmax)
ax3.set_title('Initial Distribution (Consistent Shepard Corrected Data)')
plt.colorbar(sc3, ax=ax3, label='Density')
ax3.set_xlabel('X')
ax3.set_ylabel('Y')
ax3.set_zlabel('Z')

plt.tight_layout()
plt.savefig('initial_distribution_comparison.png', dpi=300)
