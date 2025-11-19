#include "VelocityVerlet.hpp"

void velocity_verlet_step(std::vector<Particle>& particles, double dt, double dim, bool use_shepard,bool use_tensor_correction, double k, double gamma, bool use_wendland) {
    if(use_tensor_correction) {
        tensor_correction(particles, dim, use_wendland);
    }
    if(use_shepard) {
        shepard_correction(particles, dim, use_wendland);
    }

    compute_density(particles, dim, use_shepard, use_wendland);
    compute_pressure(particles, k, gamma);
    compute_sound_speed(particles, gamma);
    compute_acceleration(particles, dim, use_shepard, use_tensor_correction);
   

    // First half-step: update positions and half-step velocities
    for(auto& p : particles) {
        p.velocity = p.velocity + p.acceleration * (0.5 * dt);
        p.position = p.position + p.velocity * dt;
    }
    if(use_shepard) {
        shepard_correction(particles, dim, use_wendland);
    }   
    if(use_tensor_correction) {
        tensor_correction(particles, dim, use_wendland);
    }
    compute_density(particles, dim, use_shepard, use_wendland);
    compute_pressure(particles, k, gamma);
    compute_sound_speed(particles, gamma);
    compute_acceleration(particles, dim, use_shepard, use_tensor_correction);
    

    // Second half-step: update velocities with new accelerations
    for(auto& p : particles) {
        p.velocity = p.velocity + p.acceleration * (0.5 * dt);
    }

}