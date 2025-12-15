#include "VelocityVerlet.hpp"
void sph_leapfrog_step(
    std::vector<Particle>& particles,
    double dt,
    double dim,
    double K_0,
    double K_0_deriv,
    double G,
    double alpha_visc,
    double beta_visc,
    double epsilon_visc,
    bool use_tensor_correction,
    bool use_shepard,
    bool integrate_density)
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
     
    if(integrate_density){
        // Leapfrod density update
        density_continuity(particles, dim, use_tensor_correction);
        for(auto& p : particles){
            p.density = p.density + p.drho_dt * dt;
        }
    }
    else {
        // Recompute density from positions
        compute_density(particles, dim, use_shepard);
    }
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
    compute_acceleration(particles, dim, use_tensor_correction, alpha_visc, beta_visc, epsilon_visc);

    // -------------------------
    // 7. Full-step Kick: v_half -> v(t+dt)
    // -------------------------
    for(auto& p : particles){
        p.velocity = p.velocity + p.acceleration * 0.5 * dt;
    }
}

