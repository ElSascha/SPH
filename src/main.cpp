#include <iostream>
#include <vector>
#include "../lib/VelocityVerlet.hpp"
std::vector<Particle> load_initial_particles(char* filename){
        return std::vector<Particle>();
    }
double max_sound_speed(const std::vector<Particle>& particles){
    double max_cs = 0.0;
    for(const auto& p : particles){
        if(p.speed_of_sound > max_cs){
            max_cs = p.speed_of_sound;
        }
    }
    return max_cs;
}
int main(int argc, char* argv[]){
    std::cout << "SPH Simulation Started" << std::endl;
    // Initialize particles
    std::vector<Particle> particles = load_initial_particles(argv[1]);
    // ... (initialize your particles here)

    
    // Simulation parameters
    double CFL = 0.3;
    double steps = 100.0;
    double h = 1.0;
    double k = 1.0;
    double gamma = 1.4;
    double dt = CFL * h / max_sound_speed(particles);
    bool use_shepard = false;
    bool use_tensor_correction = false;

    // Main simulation loop
    for(int i = 0; i < steps; ++i) {
        velocity_verlet_step(particles, dt, h, 3, use_shepard, use_tensor_correction, k, gamma);
        dt = CFL * h / max_sound_speed(particles); // CFL condition update
    }

    std::cout << "SPH Simulation Ended" << std::endl;
    return 0;
}