#!/bin/bash

# Build and run script for SPH Simulation
# Usage: ./run.sh [flags for the executable]

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo "Building SPH Simulation..."

# Run make all
make all

# Check if build was successful
if [ $? -eq 0 ]; then
    echo -e "${GREEN}Build successful!${NC}"
    echo "Running simulation with flags: $@"
    echo "----------------------------------------"
    
    # Run the executable with all passed arguments
    ./build/SPH_Simulation "$@"
    
    # Store exit code
    EXIT_CODE=$?
    
    echo "----------------------------------------"
    if [ $EXIT_CODE -eq 0 ]; then
        echo -e "${GREEN}Simulation completed successfully!${NC}"
    else
        echo -e "${RED}Simulation exited with code: $EXIT_CODE${NC}"
    fi
    
    exit $EXIT_CODE
else
    echo -e "${RED}Build failed!${NC}"
    exit 1
fi
