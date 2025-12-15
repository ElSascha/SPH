#!/bin/bash
# Start an interactive allocation with srun. Options use backslash continuations
# so the command is easy to read and remains a single logical command.
srun --export=ALL \
     --time=14-00:00:00 \
     --mem-per-cpu=64G \
     --job-name=Sascha_SPH_Simulation \
     --partition=compute \
     --pty /bin/bash