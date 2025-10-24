#include "Kernel.hpp"

double W(double r, double h, int dim){
    double sigma = 0.0;
    if(dim == 1){
        sigma = (4.0/3.0) / h;
    }
    else if(dim == 2){
        sigma = (40.0/(7.0* M_PI)) / (h * h);
    }
    else if(dim == 3){
        sigma = (8.0 / M_PI) / (h*h*h);
    }
    else{
        throw std::invalid_argument("Dimension must be 1, 2, or 3.");
    }
    double q = r/h;
    double result = 0.0;
    if(0 <= q && q < (1.0/2.0)){
        result = sigma * ((6.0 * q * q * q) - (6.0 * q * q +1));
    }
    else if((1.0/2.0) <= q && q <= 1){
        result = sigma * (2 * pow(1-q,3));
    }
    else if(q > 1){
        result = 0;
    }
    else{
        throw std::invalid_argument("Wrong kernel");
    }
    return result;
}

double derivW(double r, double h, int dim){
      double sigma = 0.0;
    if(dim == 1){
        sigma = 6*(4.0/3.0) / (h*h);
    }
    else if(dim == 2){
        sigma = 6*(40.0/(7.0* M_PI)) / (h * h* h);
    }
    else if(dim == 3){
        sigma = 6*(8.0 / M_PI) / (h*h*h*h);
    }
    else{
        throw std::invalid_argument("Dimension must be 1, 2, or 3.");
    }
    double q = r / h;
    double result = 0.0;
     if(0 <= q && q < (1.0/2.0)){
        result = sigma*((3*q*q) - (2*q));
    }
    else if((1.0/2.0) <= q && q <= 1){
        result = -1.0*sigma*pow(1-q,2);
    }
    else if(q > 1){
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
    return (r/r_norm)*derivW(r_norm,h,dim);
}

void shepard_correction(std::vector<Particle>& particles, double h, int dim){
    if(dim == 3){
        for(auto& pi : particles){
            double shepard_sum = 0.0;
            for(auto& pj : particles){
                Vector r = pi.position - pj.position;
                double r_norm = r.norm();
                shepard_sum += (pj.mass / pj.density) * W(r_norm, h, dim);
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
    
void tensor_correction(std::vector<Particle>& particles, double h, int dim){
    if(dim == 3){
        for(auto& pi : particles){
            Matrix3x3 L_i;
            for(auto& pj : particles){
                Vector r = pi.position - pj.position;
                L_i = L_i + ((r * pj.mass/pj.density).outer(gradW(pi.position, pj.position, h, dim)));
            }
            try{
                pi.correction_tensor = L_i.invert();
            } catch (const std::runtime_error& e){
                pi.correction_tensor = Matrix3x3(); // set to zero matrix if not invertible
            }
        }
    }
    else{
        throw std::invalid_argument("Tensor correction only implemented for 3D.");
    }
}