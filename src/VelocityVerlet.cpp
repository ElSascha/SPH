#include "VelocityVerlet.hpp"
void sph_leapfrog_step(
    std::vector<Particle>& particles,
    double dt,
    double dim,
    double K_0,
    double K_0_deriv,
    double G,
    bool use_tensor_correction,
    bool use_shepard,
    bool use_consistent_shepard)
{
    // -------------------------
    // 1. Half-step Kick (leapfrog): v -> v_half, then x -> x + v_half*dt
    // -------------------------
    for(auto& p : particles){
        // v_half
        p.velocity = p.velocity + p.acceleration * 0.5 * dt;
        // x(t+dt) using v_half
        p.position = p.position + p.velocity * dt;
    }
    
    // -------------------------
    // 4. Update fluid state
    // -------------------------
    compute_density(particles, dim, use_shepard||use_consistent_shepard);
    if(use_tensor_correction){
        tensor_correction(particles, dim);
    }
    compute_pressure(particles, K_0, K_0_deriv);
    compute_sound_speed(particles, K_0, K_0_deriv);

    // -------------------------
    // 5. Update Stress (Hypoelastic) using current state
    //    Stress must be updated before forces that depend on it.
    // -------------------------
    compute_stress_rate(particles, G, dim, use_tensor_correction);
    for(auto& p : particles){
        p.stress = p.stress + p.dstress_dt * dt;
    }

    // -------------------------
    // 6. Compute accelerations a(t+dt)
    // -------------------------
    compute_acceleration(particles, dim, use_tensor_correction, 1.0, 2.0, 0.01);

    // -------------------------
    // 7. Full-step Kick: v_half -> v(t+dt)
    // -------------------------
    for(auto& p : particles){
        p.velocity = p.velocity + p.acceleration * 0.5 * dt;
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