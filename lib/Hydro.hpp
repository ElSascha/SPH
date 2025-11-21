#include "Kernel.hpp"
#include "Data_structs.hpp"
using namespace Data_structs;

void compute_density(std::vector<Particle>& particles, double dim, bool use_shepard, bool use_wendland);
void density_continuity(std::vector<Particle>& particles, double dim, bool use_tensor_correction, bool use_wendland);
void compute_pressure(std::vector<Particle>& particles, double k, double gamma);
void compute_acceleration(std::vector<Particle> &particles, double dim, bool use_tensor_correction, bool use_wendland);
void compute_sound_speed(std::vector<Particle>& particles, double gamma);

