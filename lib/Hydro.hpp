#include "Kernel.hpp"
#include "Data_structs.cpp"
using namespace Data_structs;

void compute_density(std::vector<Particle>& particles, double h, double dim, bool use_shepard);
void compute_pressure(std::vector<Particle>& particles, double k, double gamma);
void compute_acceleration(std::vector<Particle> &particles, double h, double dim, bool use_shepard, bool use_tensor_correction);
void compute_gravity(std::vector<Particle>& particles, double g);