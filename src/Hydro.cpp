#include "Hydro.hpp"

void compute_density(std::vector<Particle>& particles, double dim, bool use_shepard) {
    
    auto weight = [&](const Particle& pi, double r){
            if(use_wendland){
                return use_shepard ?  pi.shepard * wendland_W(r, pi.smoothing_length, dim) :  wendland_W(r, pi.smoothing_length, dim);
            } else {
                return use_shepard ?  pi.shepard * W(r, pi.smoothing_length, dim) :  W(r, pi.smoothing_length, dim);
            }
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

void density_continuity(std::vector<Particle>& particles, double dim, bool use_tensor_correction) {
    auto weight_grad = [&](const Particles& pi, const Particle& pj){
            return use_tensor_correction ? pi.correction_tensor * gradW(pi.position, pj.position, pi.smoothing_length, dim) : gradW(pi.position, pj.position, pi.smoothing_length, dim);
        };
    for(auto& pi : particles){
        pi.drho_dt = 0.0;
        for(auto& pj: particles){
            pi.drho_dt += pj.mass/pj.density * (pi.velocity - pj.velocity).dot(weight_grad(pi, pj));
        }
        pi.drho_dt *= pi.density;
    }
}

// Equation of State: Euler equation with lambda = 2
void compute_acceleration(std::vector<Particle>& particles, double dim, bool use_shepard, bool use_tensor_correction) {
        auto weight_grad = [&](const Particle& pi, const Particle& pj){
        if(use_tensor_correction == true){
            return pi.correction_tensor * gradW(pi.position, pj.position, pi.smoothing_length, dim);
        }
        else if(use_shepard == true){
            return gradW(pi.position, pj.position, pi.smoothing_length, dim) * pi.shepard;
        }
        else{
            return gradW(pi.position, pj.position, pi.smoothing_length, dim);
        }
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

// for sound speed calculation needed in time step calculation for CFL condition
void compute_sound_speed(std::vector<Particle>& particles, double gamma) {
    for(auto& pi : particles){
        pi.sound_speed = std::sqrt(gamma * pi.pressure / pi.density);
    }
}



