#!/bin/bash
#SBATCH --job-name=SPH_Simulation
#SBATCH --output=output.log
#SBATCH --error=error.log
#SBATCH --time=01:00:00
#SBATCH --partition compute
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G

# Load necessary modules
module load lib/hdf5/1.12-gnu-11.4

# Get HighFive if it doesn't exist
if [ ! -d "lib/HighFive" ]; then
    echo "Cloning HighFive library..."
    git clone https://github.com/BlueBrain/HighFive.git lib/HighFive
fi

# Build the project
make clean
make -j16
# Run the simulation
srun  build/bin/SPH_Simulation -f input/cube_distribution.h5 -N 100 
