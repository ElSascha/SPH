#include "Hydro.hpp"


void compute_density(std::vector<Particle>& particles, double dim, bool use_shepard) {
    for(auto& pi : particles){
        double density_sum = 0.0;
        const double h = pi.smoothing_length;
        const double cutoff = 2.0 * h;
        const double cutoff_sq = cutoff * cutoff;
        
        for(auto& pj : particles){
            Vector r_vec = pi.position - pj.position;
            double r_sq = r_vec.x*r_vec.x + r_vec.y*r_vec.y + r_vec.z*r_vec.z;
            
            // Early exit: skip if outside kernel support (but NOT self-contribution at r=0!)
            if(r_sq > cutoff_sq) continue;  // Changed from >= to > to include boundary
            
            double r = std::sqrt(r_sq);
            double w = W(r, h, dim);
            density_sum += pj.mass * w;
        }
        
        if(use_shepard) {
            pi.density = density_sum * pi.shepard;
        } else {
            pi.density = density_sum;
        }
    }
}

void density_continuity(std::vector<Particle>& particles, double /*dim*/) {
    for(auto& pi : particles){
        double drho_dt = 0.0;
        const double h = pi.smoothing_length;
        const double cutoff_sq = 4.0 * h * h;  // (2h)^2
        
        for(auto& pj : particles){
            if(&pi == &pj) continue;
            
            Vector r_vec = pi.position - pj.position;
            double r_sq = r_vec.x*r_vec.x + r_vec.y*r_vec.y + r_vec.z*r_vec.z;
            if(r_sq > cutoff_sq) continue;
            
            double r = std::sqrt(r_sq);
            double W_val;
            Vector gradW_ij;
            W_and_gradW(r, h, r_vec, W_val, gradW_ij);
            
            Vector vji = pj.velocity - pi.velocity;
            drho_dt += pj.mass * vji.dot(gradW_ij);
        }
        pi.drho_dt = -drho_dt;
    }
}

void compute_acceleration(std::vector<Particle>& particles, double /*dim*/, bool /*use_tensor_correction*/,
                          double alpha_visc, double beta_visc, double epsilon_visc)
{
    size_t N = particles.size();

    // Reset all accelerations to zero
    for(auto& p : particles)
        p.acceleration = Vector(0.0, 0.0, 0.0);

    for(size_t i = 0; i < N; ++i){
        const double h_i = particles[i].smoothing_length;
        const double h_j_max = h_i * 1.5;  // Assume smoothing lengths don't vary too much
        const double cutoff = 2.0 * std::max(h_i, h_j_max);
        const double cutoff_sq = cutoff * cutoff;
        
        // Cache frequently accessed values for particle i
        const Vector& pos_i = particles[i].position;
        const double rho_i = particles[i].density;
        const double rho_i_sq = rho_i * rho_i;
        const double P_i_over_rho_sq = particles[i].pressure / rho_i_sq;
        const double m_i = particles[i].mass;
        
        for(size_t j = i + 1; j < N; ++j){
            // Compute distance squared first (cheap)
            Vector r_vec = pos_i - particles[j].position;
            double r_sq = r_vec.x*r_vec.x + r_vec.y*r_vec.y + r_vec.z*r_vec.z;
            
            // Early exit: skip if outside kernel support
            if(r_sq > cutoff_sq) continue;
            
            double r = std::sqrt(r_sq);
            const double h_j = particles[j].smoothing_length;
            
            // Compute symmetrized gradient using W_and_gradW
            double W_i, W_j;
            Vector gradW_i, gradW_j;
            W_and_gradW(r, h_i, r_vec, W_i, gradW_i);
            W_and_gradW(r, h_j, r_vec, W_j, gradW_j);
            Vector gradW_ij = (gradW_i + gradW_j) * 0.5;

            // ---- Pressure + Artificial viscosity ----
            double Pi_ij = artificial_viscosity(particles[i], particles[j], alpha_visc, beta_visc, epsilon_visc);

            // Cache particle j values
            const double rho_j = particles[j].density;
            const double rho_j_sq = rho_j * rho_j;
            const double m_j = particles[j].mass;
            
            double pressure_coeff = P_i_over_rho_sq + particles[j].pressure / rho_j_sq;

            // Standard SPH momentum equation (pairwise symmetric):
            // a_i += -m_j * (P_i/rho_i^2 + P_j/rho_j^2 + Pi_ij) * gradW_ij
            // a_j += -m_i * (P_i/rho_i^2 + P_j/rho_j^2 + Pi_ij) * gradW_ji
            // Since gradW_ji = -gradW_ij:
            // a_j -= -m_i * (P_i/rho_i^2 + P_j/rho_j^2 + Pi_ij) * gradW_ij
            Vector f_ij = -(pressure_coeff + Pi_ij) * gradW_ij;

            particles[i].acceleration += f_ij * m_j;
            particles[j].acceleration -= f_ij * m_i;

            // ---- Solid stress divergence ----
            const Matrix3x3& stress_i = particles[i].stress;
            const Matrix3x3& stress_j = particles[j].stress;
            
            for(int alpha = 0; alpha < 3; ++alpha){
                double f_stress_alpha_i = 0.0;
                double f_stress_alpha_j = 0.0;
                for(int beta = 0; beta < 3; ++beta){
                    double Sij = stress_i.m[alpha][beta] / rho_i_sq +
                                 stress_j.m[alpha][beta] / rho_j_sq;

                    double gradW_beta = gradW_ij[beta];
                    f_stress_alpha_i += Sij * gradW_beta * m_j;
                    f_stress_alpha_j += Sij * gradW_beta * m_i;
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
    const double K_ratio = K_0 / K_0_deriv;
    for(auto& p: particles){
        double rho_ratio = p.density / p.rho_0;
        // Use explicit multiplication for small integer exponents
        double rho_pow = rho_ratio;
        for(int k = 1; k < static_cast<int>(K_0_deriv); ++k) {
            rho_pow *= rho_ratio;
        }
        p.pressure = K_ratio * (rho_pow - 1.0);
    }
}

void compute_sound_speed(std::vector<Particle>& particles, double K_0, double K_0_deriv)
{
    const double exp_minus_1 = K_0_deriv - 1.0;
    for(auto& p: particles){
        double rho_ratio = p.density / p.rho_0;
        // Use explicit multiplication for integer exponent (K_0_deriv - 1 = 3 for basalt)
        double rho_pow = 1.0;
        for(int k = 0; k < static_cast<int>(exp_minus_1); ++k) {
            rho_pow *= rho_ratio;
        }
        p.sound_speed = std::sqrt(K_0 * rho_pow / p.density);
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

Matrix3x3 velocity_gradient_tensor(Particle current, std::vector<Particle>& neighbors, double /*dim*/, bool use_tensor_correction){
    Matrix3x3 grad_v = Matrix3x3::zero();
    const double h = current.smoothing_length;
    const double cutoff_sq = 4.0 * h * h;
    const Vector& pos_i = current.position;
    
    if(use_tensor_correction){
        //Difference formula
        for(auto& p : neighbors){
            if(&current == &p) continue;
            
            Vector r_vec = pos_i - p.position;
            double r_sq = r_vec.x*r_vec.x + r_vec.y*r_vec.y + r_vec.z*r_vec.z;
            if(r_sq > cutoff_sq) continue;
            
            double r = std::sqrt(r_sq);
            double W_val;
            Vector gradW_ij;
            W_and_gradW(r, h, r_vec, W_val, gradW_ij);
            
            Vector vij = p.velocity - current.velocity;
            double vol_j = p.mass / p.density;
            Vector grad_psi_j_xi = current.correction_tensor * gradW_ij * vol_j;
            grad_v = grad_v + vij.outer(grad_psi_j_xi);
        }
       
    } else {
        for(auto& p : neighbors){
            if(&current == &p) continue;
            
            Vector r_vec = pos_i - p.position;
            double r_sq = r_vec.x*r_vec.x + r_vec.y*r_vec.y + r_vec.z*r_vec.z;
            if(r_sq > cutoff_sq) continue;
            
            double r = std::sqrt(r_sq);
            double W_val;
            Vector gradW_ij;
            W_and_gradW(r, h, r_vec, W_val, gradW_ij);
            
            Vector vij = p.velocity - current.velocity;
            double vol_j = p.mass / p.density;
            grad_v = grad_v + vij.outer(gradW_ij) * vol_j;
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
