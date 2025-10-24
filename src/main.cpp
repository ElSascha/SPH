#include <iostream>
#include <vector>
#include "../lib/VelocityVerlet.hpp"

int main(){
    std::cout << "SPH Simulation Started" << std::endl;

    // Initialize particles
    std::vector<Particle> particles;
    // ... (initialize your particles here)

    // Simulation parameters
    double dt = 0.01;
    double h = 1.0;
    double k = 1.0;
    double gamma = 1.4;
    bool use_shepard = true;

    // Main simulation loop
    for(int i = 0; i < 100; ++i) {
        velocity_verlet_step(particles, dt, h, 3, use_shepard, k, gamma);
    }

    std::cout << "SPH Simulation Ended" << std::endl;
    return 0;
}