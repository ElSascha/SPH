#include "Hydro.hpp"

void compute_density(std::vector<Particle>& particles, double h, double dim, bool use_shepard) {
    for(auto& pi : particles){
        double density_sum = 0.0;
        pi.density = 0.0;
        for(auto& pj : particles){
            double r_norm = (pi.position - pj.position).norm();
            if(use_shepard){
                density_sum += pj.mass *  pi.shepard * W(r_norm, h, dim);
            }
            else{
                density_sum += pj.mass * W(r_norm, h, dim);
            }
        }
        pi.density = density_sum;
    }
}
// Equation of State: Euler equation with lambda = 2
void compute_acceleration(std::vector<Particle>& particles, double h, double dim, bool use_shepard) {
 
}
