#include "VelocityVerlet.hpp"
void sph_leapfrog_step(
    std::vector<Particle>& particles,
    double dt,
    double dim,
    double K_0,
    double K_0_deriv,
    double rho_0,
    double G,
    bool use_tensor_correction,
    bool use_shepard)
{
    // -------------------------
    // 1. Half-step Kick: v(t) -> v(t + dt/2)
    // -------------------------
    std::vector<Vector> old_velocity = std::vector<Vector>(particles.size()); 
    for(auto& p : particles){
        old_velocity[&p - &particles[0]] = p.velocity; // store old velocity
        p.position = p.position + p.velocity * 0.5 * dt;
        p.velocity = p.velocity + p.acceleration * 0.5 * dt;
    }
    
    // -------------------------
    // 4. Update fluid state
    // -------------------------
    compute_density(particles, dim, use_shepard);
    if(use_tensor_correction){
        tensor_correction(particles, dim);
    }
    compute_pressure(particles, K_0, K_0_deriv, rho_0);
    compute_sound_speed(particles, K_0, K_0_deriv, rho_0);

    // Update Stress (Hypoelastic)
    compute_stress_rate(particles, G, dim, use_tensor_correction);
    for(auto& p : particles){
        p.stress = p.stress + p.dstress_dt * dt;
    }

    // -------------------------
    // 5. Compute accelerations a(t+dt)
    // -------------------------
    compute_acceleration(particles, dim, use_tensor_correction, 1.0, 2.0, 0.01);

    // -------------------------
    // 6. Full-step Kick: v(t + dt/2) -> v(t + dt)
    // -------------------------
    for(auto& p : particles){
        p.velocity = old_velocity[&p - &particles[0]] + p.acceleration * dt;
        p.position =  p.position + p.velocity * 0.5 * dt;
    }
}

void integrate_density(std::vector<Particle>& particles, double dt, double dim)
{
    // ------ 1) Compute drho_dt(t) ------ 
    density_continuity(particles, dim);

    std::vector<double> drho_dt_old(particles.size());
    std::vector<double> rho_old(particles.size());

    for (size_t i = 0; i < particles.size(); ++i) {
        rho_old[i] = particles[i].density;
        drho_dt_old[i] = particles[i].drho_dt;
        particles[i].density = rho_old[i] + dt * drho_dt_old[i];   // predictor
    }

    // ------ 2) drho_dt at predicted state ------
    density_continuity(particles, dim);

    // ------ 3) Corrector step ------
    for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].density =
            rho_old[i] + 0.5 * dt * (drho_dt_old[i] + particles[i].drho_dt);
    }

}