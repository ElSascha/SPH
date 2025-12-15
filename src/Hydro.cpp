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

void density_continuity(std::vector<Particle>& particles, double /*dim*/, bool use_tensor_correction){ 
    for(auto& pi : particles){
        double nabla_dot_v = grad_dot_v(pi, particles, 3, use_tensor_correction);
        pi.drho_dt = -pi.density * nabla_dot_v;   
    }
}

void compute_acceleration(std::vector<Particle>& particles, double /*dim*/, bool use_tensor_correction,
                          double alpha_visc, double beta_visc, double epsilon_visc)
{
    size_t N = particles.size();

    // Reset all accelerations to zero
    for(auto& p : particles)
        p.acceleration = Vector(0.0, 0.0, 0.0);

    for(size_t i = 0; i < N; ++i){
        const double h_i = particles[i].smoothing_length;
        // Optimization: Pre-calculate constants for particle i
        const Vector& pos_i = particles[i].position;
        const double rho_i = particles[i].density;
        const double rho_i_sq = rho_i * rho_i;
        const double P_i_term = particles[i].pressure / rho_i_sq;
        const Matrix3x3 stress_i_term = particles[i].stress * (1.0 / rho_i_sq);
        const double m_i = particles[i].mass;

        // Neighbor search parameters
        const double h_j_max_est = h_i * 1.5; // Estimation for cutoff
        const double cutoff = 2.0 * std::max(h_i, h_j_max_est);
        const double cutoff_sq = cutoff * cutoff;
        
        for(size_t j = i + 1; j < N; ++j){
            Vector r_vec = pos_i - particles[j].position;
            double r_sq = r_vec.x*r_vec.x + r_vec.y*r_vec.y + r_vec.z*r_vec.z;
            
            if(r_sq > cutoff_sq) continue;
            
            double r = std::sqrt(r_sq);
            const double h_j = particles[j].smoothing_length;
            
            // 1. Calculate Raw Gradient (Symmetrized Kernel)
            double W_i, W_j;
            Vector gradW_raw_i, gradW_raw_j;
            W_and_gradW(r, h_i, r_vec, W_i, gradW_raw_i);
            W_and_gradW(r, h_j, r_vec, W_j, gradW_raw_j);
            
            // Standard SPH gradient of W_ij w.r.t x_i
            Vector gradW_raw = (gradW_raw_i + gradW_raw_j) * 0.5;

            // 2. Apply Tensor Correction (KGC)
            // If enabled, applying L_i to i's terms and L_j to j's terms
            Vector gradW_i_term, gradW_j_term;
            if (use_tensor_correction) {
                gradW_i_term = particles[i].correction_tensor * gradW_raw;
                gradW_j_term = particles[j].correction_tensor * gradW_raw;
            } else {
                gradW_i_term = gradW_raw;
                gradW_j_term = gradW_raw;
            }

            // 3. Pressure & Viscosity Force (Symmetric)
            const double rho_j = particles[j].density;
            const double rho_j_sq = rho_j * rho_j;
            const double m_j = particles[j].mass;
            
            double Pi_ij = artificial_viscosity(particles[i], particles[j], alpha_visc, beta_visc, epsilon_visc);
            
            // Distribute viscosity symmetrically (0.5 to each term)
            double term_scalar_i = P_i_term + (0.5 * Pi_ij);
            double term_scalar_j = (particles[j].pressure / rho_j_sq) + (0.5 * Pi_ij);

            // Compute separate vectors to handle the different corrections
            Vector force_vec_i = gradW_i_term * term_scalar_i;
            Vector force_vec_j = gradW_j_term * term_scalar_j;

            // Total Pressure Force F_ij (Notice negative sign for pressure gradient)
            Vector f_pressure = (force_vec_i + force_vec_j) * (-1.0);

            particles[i].acceleration += f_pressure * m_j;
            particles[j].acceleration -= f_pressure * m_i;

            // 4. Solid Stress Force (Symmetric)
            const Matrix3x3 stress_j_term = particles[j].stress * (1.0 / rho_j_sq);
            
            for(int alpha = 0; alpha < 3; ++alpha){
                double f_stress_scalar = 0.0;
                
                // Term from Particle i (uses L_i corrected gradient)
                for(int beta = 0; beta < 3; ++beta){
                    f_stress_scalar += stress_i_term.m[alpha][beta] * gradW_i_term[beta];
                }

                // Term from Particle j (uses L_j corrected gradient)
                for(int beta = 0; beta < 3; ++beta){
                    f_stress_scalar += stress_j_term.m[alpha][beta] * gradW_j_term[beta];
                }
                
                particles[i].acceleration[alpha] += f_stress_scalar * m_j;
                particles[j].acceleration[alpha] -= f_stress_scalar * m_i; 
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


double grad_dot_v(Particle current, std::vector<Particle>& neighbors, double /*dim*/, bool use_tensor_correction){
    const double h = current.smoothing_length;
    const double cutoff_sq = 4.0 * h * h;
    const Vector& pos_i = current.position;
    
    double result = 0.0;
    
    if(use_tensor_correction){
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
            result += vij.dot(grad_psi_j_xi);
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
            result += vij.dot(gradW_ij) * vol_j;
        }
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
