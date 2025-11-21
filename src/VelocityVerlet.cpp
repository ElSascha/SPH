#include "VelocityVerlet.hpp"
void sph_leapfrog_step(
    std::vector<Particle>& particles,
    double dt,
    double dim,
    double k,
    double gamma,
    bool use_tensor_correction,
    bool use_wendland)
{
    shepard_correction(particles, dim, use_wendland);
    consistent_shepard_interpolation(particles, dim, use_wendland);
    tensor_correction(particles, dim, use_wendland);
    // -------------------------
    // 0. Kick (t -> t + dt/2)
    // -------------------------
    for (auto& p : particles)
        p.velocity +=  p.acceleration * 0.5 * dt;

    // -------------------------
    // 1. Drift positions
    // -------------------------
    for (auto& p : particles)
        p.position += p.velocity * dt;


    // -------------------------
    // 2. Update density
    //    (continuity equation)
    // -------------------------
    density_update(particles, dt, dim, use_wendland, use_tensor_correction);

    // -------------------------
    // 3. Recompute fluid state
    // -------------------------
    compute_pressure(particles, k, gamma);
    compute_sound_speed(particles, gamma);

    // -------------------------
    // 6. Compute accelerations a(t + dt/2)
    // -------------------------
    compute_acceleration(particles, dim, use_tensor_correction, use_wendland);

    // -------------------------
    // 7. Kick (t + dt/2 -> t + dt)
    // -------------------------
    for (auto& p : particles)
        p.velocity += p.acceleration * 0.5 * dt;

    
}

void density_update(
    std::vector<Particle>& particles,
    double dt,
    double dim,
    bool use_wendland,
    bool use_tensor_correction
    )
{
    // predictor
    for (auto& p : particles)
        p.drho_dt = 0.0;

    density_continuity(particles, dim, use_tensor_correction, use_wendland);

    for (auto& p : particles)
        p.rho_pred = p.density + 0.5 * dt * p.drho_dt;

    // corrector
    // Temporarily update density to predicted density
    for (auto& p : particles)
        p.density = p.rho_pred;

    for (auto& p : particles)
        p.drho_dt = 0.0;

    density_continuity(particles, dim, use_tensor_correction, use_wendland);

    // Restore original density before the final update
    for (auto& p : particles)
        p.density = p.rho_pred - 0.5 * dt * p.drho_dt;

    for (auto& p : particles)
        p.density += dt * p.drho_dt;
    
}