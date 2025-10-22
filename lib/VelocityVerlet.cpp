#include "VelocityVerlet.hpp"

void velocity_verlet_step(std::vector<Particle>& particles, double dt, double h, double dim, bool use_shepard, double k, double gamma) {
    compute_acceleration(particles, h, dim, use_shepard);
    compute_density(particles, h, dim, use_shepard);
    compute_pressure(particles, k, gamma);

    // First half-step: update positions and half-step velocities
    for(auto& p : particles) {
        p.velocity = p.velocity + p.acceleration * (0.5 * dt);
        p.position = p.position + p.velocity * dt;
    }
    compute_acceleration(particles, h, dim, use_shepard);
    compute_density(particles, h, dim, use_shepard);
    compute_pressure(particles, k, gamma);


    // Second half-step: update velocities with new accelerations
    for(auto& p : particles) {
        p.velocity = p.velocity + p.acceleration * (0.5 * dt);
    }

}