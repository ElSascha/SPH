#include "Hydro.hpp"


void compute_density(std::vector<Particle>& particles, double dim, bool use_shepard) {
    auto weight = [&](const Particle& pi, double r){
        return use_shepard ?  pi.shepard * W(r, pi.smoothing_length, dim) :  W(r, pi.smoothing_length, dim);
    };
    for(auto& pi : particles){
        double density_sum = 0.0;
        for(auto& pj : particles){
            double r_norm = (pi.position - pj.position).norm();
            density_sum += pj.mass * weight(pi, r_norm);
        }
        pi.density = density_sum;
    }
}

void density_continuity(std::vector<Particle>& particles, double dim) {
    for(auto& pi : particles){
        double drho_dt = 0.0;
        for(auto& pj : particles){
            if(&pi == &pj) continue;
            Vector vji = pj.velocity - pi.velocity;
            Vector gradW_ij = gradW(pi.position, pj.position, pi.smoothing_length, dim);
            drho_dt += pj.mass * vji.dot(gradW_ij);
        }
        pi.drho_dt = -drho_dt;
    }
}

void compute_acceleration(std::vector<Particle>& particles, double dim, bool use_tensor_correction,
                          double alpha_visc, double beta_visc, double epsilon_visc)
{
    size_t N = particles.size();

    // Zuerst alle Beschleunigungen auf null setzen
    for(auto& p : particles)
        p.acceleration = Vector(0.0, 0.0, 0.0);

    for(size_t i = 0; i < N; ++i){
        for(size_t j = i + 1; j < N; ++j){  // nur j>i für Paarweise
            if(i == j) continue;
            
            // Symmetrized kernel gradient for angular momentum conservation
            // Use average of gradients computed with both smoothing lengths
            Vector gradW_i = gradW(
                particles[i].position,
                particles[j].position,
                particles[i].smoothing_length,
                dim
            );
            Vector gradW_j = gradW(
                particles[i].position,
                particles[j].position,
                particles[j].smoothing_length,
                dim
            );
            Vector gradW_ij = (gradW_i + gradW_j) * 0.5;

            // ---- Pressure + Artificial viscosity ----
            double Pi_ij = artificial_viscosity(particles[i], particles[j], alpha_visc, beta_visc, epsilon_visc);

            // Symmetric pressure term: P_i/rho_i^2 + P_j/rho_j^2
            double pressure_coeff = particles[i].pressure / (particles[i].density * particles[i].density) +
                                    particles[j].pressure / (particles[j].density * particles[j].density);

            // Symmetric formulation: each particle uses its own mass
            // a_i += -m_j * (P_i/rho_i^2 + P_j/rho_j^2 + Pi_ij) * gradW_ij
            // a_j -= -m_i * (P_i/rho_i^2 + P_j/rho_j^2 + Pi_ij) * gradW_ij  (note: gradW_ji = -gradW_ij)
            Vector f_pressure_i = -(pressure_coeff + Pi_ij) * gradW_ij * particles[j].mass;
            Vector f_pressure_j = -(pressure_coeff + Pi_ij) * gradW_ij * particles[i].mass;

            particles[i].acceleration += f_pressure_i;
            particles[j].acceleration -= f_pressure_j;

            // ---- Solid stress divergence ----
            // Symmetric stress coefficient
            double m_i = particles[i].mass;
            double m_j = particles[j].mass;
            
            for(int alpha = 0; alpha < 3; ++alpha){
                double f_stress_alpha_i = 0.0;
                double f_stress_alpha_j = 0.0;
                for(int beta = 0; beta < 3; ++beta){
                    double Sij = particles[i].stress.m[alpha][beta] / (particles[i].density * particles[i].density) +
                                 particles[j].stress.m[alpha][beta] / (particles[j].density * particles[j].density);

                    f_stress_alpha_i += Sij * gradW_ij[beta] * m_j;
                    f_stress_alpha_j += Sij * gradW_ij[beta] * m_i;
                }
                particles[i].acceleration[alpha] += f_stress_alpha_i;
                particles[j].acceleration[alpha] -= f_stress_alpha_j;
            }
        }
    }
}               

double artificial_viscosity(const Particle& pi, const Particle& pj, double alpha, double beta, double epsilon) {
    Vector rij = pi.position - pj.position;
    Vector vij = pi.velocity - pj.velocity;

    double rij_dot_vij = vij.dot(rij);
    if(rij_dot_vij >= 0.0) return 0.0;

    double r2 = std::max(rij.dot(rij), 1e-24);
    double h_avg = 0.5 * (pi.smoothing_length + pj.smoothing_length);

    double mu_ij = (h_avg * rij_dot_vij) / (r2 + epsilon * h_avg * h_avg);

    double c = std::max(0.5 * (pi.sound_speed + pj.sound_speed), 1e-12);
    double rho = 0.5 * (pi.density + pj.density);

    double Pi_ij = (-alpha * c * mu_ij + beta * mu_ij * mu_ij) / rho;
    return Pi_ij;
}


// Murnaghan EOS: P = K/(gamma) * [(rho/rho0)^gamma - 1]; c_s = sqrt(dP/drho)
void compute_pressure(std::vector<Particle>& particles, double K_0, double K_0_deriv)
{
    for(auto& p: particles){
        p.pressure = K_0 / K_0_deriv * (pow(p.density/p.rho_0, K_0_deriv) - 1.0);
        // Optionally clamp to prevent excessive tension (uncomment if needed):
        // p.pressure = std::max(p.pressure, -0.1 * K_0);  // Limit tension to 10% of bulk modulus
    }
}

void compute_sound_speed(std::vector<Particle>& particles, double K_0, double K_0_deriv)
{
    for(auto& p: particles){
        p.sound_speed = std::sqrt(K_0 * pow(p.density/p.rho_0,K_0_deriv - 1) / p.density);
    }
}

// Solid mechanics specific functions 
double compute_shear_modulus(double E, double nu){
    return E / (2.0 * (1.0 + nu));
}



double grad_dot_v(Particle current, std::vector<Particle>& neighbors, double dim, bool use_tensor_correction){
    if(use_tensor_correction){
        return 0.0; // placeholder
    }
    double result = 0.0;
    for(auto& p : neighbors){
        if(&current == &p) continue;
        Vector vij = p.velocity - current.velocity;
        Vector gradW_ij = gradW(current.position, p.position, current.smoothing_length, dim);
        result += p.mass / p.density * vij.dot(gradW_ij);
    }
    return result;    
}

Matrix3x3 velocity_gradient_tensor(Particle current, std::vector<Particle>& neighbors, double dim, bool use_tensor_correction){
    Matrix3x3 grad_v = Matrix3x3::zero(); // initialize zero matrix
    if(use_tensor_correction){
        // implement the following scheme: (v_j - v_i) ∇ psi_j(x_i)^T
        // where ∇ psi_j(x_i) = m_j/rho_j * L_i ∇ W_ij
        for(auto& p : neighbors){
            if(&current == &p) continue;
            Vector vij = p.velocity - current.velocity;
            Vector grad_psi_j_xi = current.correction_tensor * gradW(current.position, p.position, current.smoothing_length, dim) * (p.mass / p.density);
            grad_v =  grad_v + vij.outer(grad_psi_j_xi);
        }
    } else {
        // Standard SPH gradient without tensor correction
        for(auto& p : neighbors){
            if(&current == &p) continue;
            Vector vij = p.velocity - current.velocity;
            Vector gradW_ij = gradW(current.position, p.position, current.smoothing_length, dim);
            grad_v = grad_v + vij.outer(gradW_ij) * (p.mass / p.density);
        }
    }
    return grad_v;
}

void compute_stress_rate(std::vector<Particle>& particles, double G, double dim, bool use_tensor_correction) {
    for(auto& p : particles) {
        Matrix3x3 grad_v = velocity_gradient_tensor(p, particles, dim, use_tensor_correction);
        Matrix3x3 strain_rate = (grad_v + grad_v.transpose()) * 0.5;
        Matrix3x3 spin = (grad_v - grad_v.transpose()) * 0.5;
        
        // Deviatoric strain rate
        double trace = (strain_rate.m[0][0] + strain_rate.m[1][1] + strain_rate.m[2][2]) / 3.0;
        Matrix3x3 dev_strain_rate = strain_rate;
        for(int i=0; i<3; ++i) dev_strain_rate.m[i][i] -= trace;

        // Jaumann rate: dsigma/dt = sigma_obj_rate + Omega * sigma - sigma * Omega
        // sigma_obj_rate = 2G * dev_epsilon_dot
        Matrix3x3 sigma = p.stress;
        Matrix3x3 term1 = spin * sigma;
        Matrix3x3 term2 = sigma * spin * (-1.0); // -sigma * Omega
        Matrix3x3 term3 = dev_strain_rate * (2.0 * G);

        p.dstress_dt = term1 + term2 + term3;
    }
}
