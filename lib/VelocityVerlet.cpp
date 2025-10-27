#include "VelocityVerlet.hpp"

void velocity_verlet_step(std::vector<Particle>& particles, double dt, double h, double dim, bool use_shepard,bool use_tensor_correction, double k, double gamma) {
    if(use_tensor_correction) {
        tensor_correction(particles, h, dim);
    }
    if(use_shepard) {
        shepard_correction(particles, h, dim);
    }

    compute_density(particles, h, dim, use_shepard);
    compute_pressure(particles, k, gamma);
    compute_sound_speed(particles, k, gamma);
    compute_acceleration(particles, h, dim, use_shepard, use_tensor_correction);
    compute_gravity(particles, 9.81); // Example gravity value
   

    // First half-step: update positions and half-step velocities
    for(auto& p : particles) {
        p.velocity = p.velocity + p.acceleration * (0.5 * dt);
        p.position = p.position + p.velocity * dt;
    }
    if(use_shepard) {
        shepard_correction(particles, h, dim);
    }   
    compute_density(particles, h, dim, use_shepard);
    compute_pressure(particles, k, gamma);
    compute_sound_speed(particles, k, gamma);
    compute_acceleration(particles, h, dim, use_shepard, use_tensor_correction);
    compute_gravity(particles, 9.81); // Example gravity value

    // Second half-step: update velocities with new accelerations
    for(auto& p : particles) {
        p.velocity = p.velocity + p.acceleration * (0.5 * dt);
    }

}