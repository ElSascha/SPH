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

void consistent_shepard_correction(std::vector<Particle>& particles, double h, int dim){
    for(auto& pi : particles){
        
    }
}
