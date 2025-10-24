#include "Kernel.hpp"
#include "Data_structs.hpp"
using namespace Data_structs;

void compute_density(std::vector<Particle>& particles, double h, double dim, bool use_shepard);
void compute_pressure(std::vector<Particle>& particles, double k, double gamma);
void compute_acceleration(std::vector<Particle> &particles, double h, double dim, bool use_shepard, bool use_tensor_correction);
void compute_gravity(std::vector<Particle>& particles, double g);
void compute_vector_field(std::vector<Particle>& particles, double h, double dim, bool use_shepard); // for testing purposes
void compute_tensor_field(std::vector<Particle>& particles, double h, double dim, bool use_tensor_correction); // for testing purposes