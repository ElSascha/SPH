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
    // Tensor nicht hier verwenden. Tensor gehört NACH dieser Integration (siehe main loop).
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

    // 2. Berechne Beschleunigungen
    for(size_t i=0; i<N; ++i){
        Vector acc(0.0, 0.0, 0.0);

        for(size_t j=0; j<N; ++j){
            if(i == j) continue;
            Vector gradW_ij = gradW(particles[i].position, particles[j].position, particles[i].smoothing_length, dim);

            // --- Druckterm + artificial viscosity ---
            double Pi_ij = artificial_viscosity(particles[i], particles[j], alpha_visc, beta_visc, epsilon_visc);
            double pressure_term = particles[i].pressure/(particles[i].density*particles[i].density) +
                                   particles[j].pressure/(particles[j].density*particles[j].density) +
                                   Pi_ij;
            Vector f_pressure = -particles[j].mass * pressure_term * gradW_ij;

            // --- Shearterm (Matrix-Vektor-Produkt) ---
            // Uses the pre-calculated stress tensor from the Hypoelastic update
            Vector f_shear(0.0, 0.0, 0.0);
            for(int alpha=0; alpha<3; ++alpha){
                for(int beta=0; beta<3; ++beta){
                    f_shear[alpha] += (particles[i].stress.m[alpha][beta]/(particles[i].density*particles[i].density) +
                                       particles[j].stress.m[alpha][beta]/(particles[j].density*particles[j].density)) *
                                       gradW_ij[beta] * particles[j].mass;
                }
            }

            acc += f_pressure + f_shear;
        }

        particles[i].acceleration = acc;
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
void compute_pressure(std::vector<Particle>& particles, double K_0, double K_0_deriv, double rho_0)
{
    for(auto& p: particles){
        p.pressure = K_0 / K_0_deriv *(pow(p.density/rho_0,K_0_deriv) - 1.0);
    }
}

void compute_sound_speed(std::vector<Particle>& particles, double K_0, double K_0_deriv, double rho_0)
{
    for(auto& p: particles){
        p.sound_speed = std::sqrt(K_0 * pow(p.density/rho_0,K_0_deriv - 1) / p.density);
    }
}

// Solid mechanics specific functions 
double compute_shear_modulus(double E, double nu){
    return E / (2.0 * (1.0 + nu));
}



double grad_v(Particle current, std::vector<Particle>& neighbors, double dim, bool use_tensor_correction){
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
        return grad_v; // placeholder
    }
    for(auto& p : neighbors){
        if(&current == &p) continue;
        Vector vij = p.velocity - current.velocity;
        Vector gradW_ij = gradW(current.position, p.position, current.smoothing_length, dim);
        grad_v = grad_v + vij.outer(gradW_ij) * (p.mass / p.density);
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

        // Jaumann rate: dsigma/dt = sigma * Omega^T + Omega * sigma + 2G * dev_epsilon_dot
        // Note: Omega^T = -Omega
        Matrix3x3 sigma = p.stress;
        Matrix3x3 term1 = sigma * spin.transpose();
        Matrix3x3 term2 = spin * sigma;
        Matrix3x3 term3 = dev_strain_rate * (2.0 * G);

        p.dstress_dt = term1 + term2 + term3;
    }
}

// Compute Strain
Matrix3x3 compute_strain_tensor(const Matrix3x3& grad_v){
    Matrix3x3 strain = (grad_v + grad_v.transpose()) * 0.5;
    return strain;
}

Matrix3x3 compute_deviatoric_strain(const Matrix3x3& strain){
    double trace = (strain.m[0][0] + strain.m[1][1] + strain.m[2][2]) / 3.0;
    Matrix3x3 dev_strain = strain;
    for(int i=0;i<3;++i){
        dev_strain.m[i][i] -= trace;
    }
    return dev_strain;
}

Matrix3x3 compute_shear_stress(const Matrix3x3& dev_strain, double G){
    Matrix3x3 shear_stress = dev_strain * (2.0 * G);
    return shear_stress;
}
