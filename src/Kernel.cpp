#include "Kernel.hpp"

double W(double r, double h, int dim){
    double sigma = 0.0;
    if(dim == 1){
        sigma = (2.0/3.0) / h;
    }
    else if(dim == 2){
        sigma = (10.0/(7.0* M_PI)) / (h * h);
    }
    else if(dim == 3){
        sigma = (1.0 / M_PI) / (h*h*h);
    }
    else{
        throw std::invalid_argument("Dimension must be 1, 2, or 3.");
    }
    double q = r/h;
    double result = 0.0;
    if(q >= 0 && q < 1){
        result = sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q);
    }
    else if(q >= 1 && q < 2){
        result = sigma * (0.25 * (2.0-q)*(2.0-q)*(2.0-q));
    }
    else if(q >= 2){
        result = 0.0;
    }
    else{
        throw std::invalid_argument("Wrong kernel");
    }
    return result;
}

double wendland_W(double r, double h, int dim){
    double sigma = 0.0;
    double q = r / h;
    if(dim == 3){
        sigma = 495.0 / (256.0 * M_PI * h*h*h);
    }
    else{
        throw std::invalid_argument("Wendland kernel only implemented for 3D.");
    }
    if(q < 2.0){
        double term = 1.0 - q/2.0;
        return sigma * pow(term, 6.0) * (1.0 + 6.0 * (q/2.0) + 35.0/3.0 * (q/2.0)*(q/2.0));
    }
    else{
        return 0.0;
    }
}

double wendland_W_deriv(double r, double h, int dim){
    double sigma = 0.0;
    double q = r / h;
    if(dim == 3){
        sigma = 495.0 / (256.0 * M_PI * h*h*h);
    }
    else{
        throw std::invalid_argument("Wendland kernel only implemented for 3D.");
    }
    if(q < 2.0){
        double term = 1.0 - q/2.0;
        return sigma * ( -3.0 * pow(term, 5.0) * (1.0 + 6.0 * (q/2.0) + 35.0/3.0 * (q/2.0)*(q/2.0)) + pow(term, 6.0) * (3.0 + 35.0/3.0 * q/2.0) ) / h;
    }
    else{
        return 0.0;
    }
}

double derivW(double r, double h, int dim){
    double sigma = 0.0;
    if(dim == 1){
        sigma = (2.0/3.0) / (h*h);
    }
    else if(dim == 2){
        sigma = (10.0/(7.0* M_PI)) / (h*h*h);
    }
    else if(dim == 3){
        sigma = (1.0 / M_PI) / (h*h*h*h);
    }
    else{
        throw std::invalid_argument("Dimension must be 1, 2, or 3.");
    }
    double q = r / h;
    double result = 0.0;
     if(q >= 0 && q < 1){
        result = sigma*(-3.0*q + 2.25*q*q);
    }
    else if(q >= 1 && q < 2){
        result = -0.75*sigma*(2.0-q)*(2.0-q);
    }
    else if(q >= 2){
        result = 0;
    }
    else{
        throw std::invalid_argument("Wrong kernel");
    }
    return result;


}

Vector gradW(Vector pos_i, Vector pos_j, double h, int dim){
    Vector r = pos_i - pos_j;
    double r_norm = r.norm();
    if (r_norm == 0.0) {
        return Vector(0,0,0);
    }
    return (r/r_norm)*derivW(r_norm,h,dim);
}

void shepard_correction(std::vector<Particle>& particles, int dim, bool use_wendland){
    if(dim == 3){
        if(use_wendland){
            for(auto& pi : particles){
            double shepard_sum = 0.0;
            for(auto& pj : particles){
                Vector r = pi.position - pj.position;
                double r_norm = r.norm();
                shepard_sum += (pj.mass / pj.density) * wendland_W(r_norm, pj.smoothing_length, dim);
            }
            if(shepard_sum > 1e-12){
                pi.shepard = 1.0/shepard_sum;
            }
            else{
                pi.shepard = 1.0;
            }
        } 
            return;
        }
        for(auto& pi : particles){
            double shepard_sum = 0.0;
            for(auto& pj : particles){
                Vector r = pi.position - pj.position;
                double r_norm = r.norm();
                shepard_sum += (pj.mass / pj.density) * W(r_norm, pj.smoothing_length, dim);
            }
            if(shepard_sum > 1e-12){
                pi.shepard = 1.0/shepard_sum;
            }
            else{
                pi.shepard = 1.0;
            }
        } 
    }
    else{
        throw std::invalid_argument("Shepard correction only implemented for 3D.");
    }
}


void consistent_shepard_interpolation(std::vector<Particle>& particles, int dim, bool use_wendland){
    if(dim != 3){
        throw std::invalid_argument("Shepard correction only implemented for 3D.");
    }

    const int max_iter = 50;
    const double tol = 1e-6;
    // set all shepard factors to 1.0 initially
    for(auto& p : particles) p.shepard = 1.0;
    // new and old shepard factors for convergence check
    std::vector<double> c_old(particles.size(), 1.0);
    std::vector<double> c_new(particles.size(), 1.0);
    // power method iteration
    for(int iter = 0; iter < max_iter; ++iter){
       for(size_t i = 0; i < particles.size(); ++i){
            double shepard_sum = 0.0;
            auto& pi = particles[i];
            for(size_t j = 0; j < particles.size(); ++j){
                auto& pj = particles[j];
                Vector r = pi.position - pj.position;
                double r_norm = r.norm();
                if(use_wendland){
                    shepard_sum += (pj.mass / (c_old[j] * pj.density)) * wendland_W(r_norm, pj.smoothing_length, dim);
                } else {
                    shepard_sum += (pj.mass / (c_old[j] * pj.density)) * W(r_norm, pj.smoothing_length, dim);
                }
            }
            if(shepard_sum > 1e-12){
                c_new[i] = 1.0 / shepard_sum;
            }
            else{
                c_new[i] = 1.0;
            }
       }
       // check for convergence
       double max_diff = 0.0;
       for(size_t i = 0; i < particles.size(); ++i){
           max_diff = std::max(max_diff, std::abs(c_new[i] - c_old[i]));
       }
       if(max_diff < tol){
            
           break;
       }
       c_old = c_new;
    }
    // update shepard factors
    for(size_t i = 0; i < particles.size(); ++i){
        particles[i].shepard = c_new[i];
    }
}

    
void tensor_correction(std::vector<Particle>& particles, int dim, bool use_wendland){
    if(dim == 3){
    
        for(auto& pi : particles){
            Matrix3x3 L_i;
            for(auto& pj : particles){
                Vector r = pi.position - pj.position;
                if(use_wendland){
                    L_i = L_i + ((r * pj.mass/pj.density).outer(gradW(pi.position, pj.position, pj.smoothing_length, dim)));
                } else {
                    L_i = L_i + ((r * pj.mass/pj.density).outer(gradW(pi.position, pj.position, pj.smoothing_length, dim)));
                }
            }
            try{
                pi.correction_tensor = L_i.invert();
            } catch (const std::runtime_error& e){
                pi.correction_tensor = id(); // set to identity matrix if not invertible
            }
        }
    }
    else{
        throw std::invalid_argument("Tensor correction only implemented for 3D.");
    }
}