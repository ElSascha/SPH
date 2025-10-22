#include "Hydro.hpp"

void compute_density(std::vector<Particle>& particles, double h, double dim, bool use_shepard) {
    auto weight = [&](const Particle pi, double r){
            return use_shepard ?  pi.shepard * W(r, h, dim) :  W(r, h, dim);
        };
    for(auto& pi : particles){
        double density_sum = 0.0;
        pi.density = 0.0;
        for(auto& pj : particles){
            double r_norm = (pi.position - pj.position).norm();
            density_sum += pj.mass * weight(pi, r_norm);
        }
        pi.density = density_sum;
    }
}
// Equation of State: Euler equation with lambda = 2
void compute_acceleration(std::vector<Particle>& particles, double h, double dim, bool use_shepard) {
    auto weight_grad = [&](const Particle pi, const Particle pj){
        return use_shepard? gradW(pi.position, pj.position, h, dim) * pi.shepard : gradW(pi.position, pj.position, h, dim);
    };
    for(auto& pi : particles){
        Vector acc_sum(0,0,0);
        pi.acceleration = Vector(0,0,0);
        for(auto& pj : particles){
            if(&(pi) == &(pj)) continue; // skip self-interaction
            Vector gradW_ij = weight_grad(pi, pj);
            Vector pressure_term = gradW_ij * (pi.pressure / (pi.density * pi.density) + pj.pressure / (pj.density * pj.density));
            acc_sum =  pressure_term + acc_sum;
        }
        pi.acceleration = acc_sum * (-1.0);
    }
    }
    
    void compute_pressure(std::vector<Particle>& particles, double k, double gamma) {
    for(auto& pi : particles){
        pi.pressure = k * pow(pi.density, gamma);
    }
    }

    // for surface stability test
    void compute_gravity(std::vector<Particle>& particles, double g) {
    for(auto& pi : particles){
        Vector gravity = Vector(0, -g, 0);
        pi.acceleration = pi.acceleration + gravity;
    }
}
