#pragma once
#include "Kernel.hpp"
#include "Data_structs.hpp"
using namespace Data_structs;

void compute_density(std::vector<Particle>& particles, double dim, bool use_shepard);
void density_continuity(std::vector<Particle>& particles, double dim);
void compute_acceleration(std::vector<Particle>& particles, double dim, bool use_tensor_correction,
                          double alpha_visc=1.0, double beta_visc=2.0, double epsilon_visc=0.01);
double artificial_viscosity(const Particle& pi, const Particle& pj, double alpha, double beta, double epsilon);
// Murnaghan EOS: P = K/(gamma) * [(rho/rho0)^gamma - 1]; c_s = sqrt(dP/drho)
void compute_pressure(std::vector<Particle>& particles, double K_0, double K_0_deriv);
void compute_sound_speed(std::vector<Particle>& particles, double K_0, double K_0_deriv);
// Solid mechanics specific functions 
double compute_shear_modulus(double E, double nu);
double grad_dot_v(Particle current, std::vector<Particle>& neighbors, double dim, bool use_tensor_correction);
Matrix3x3 velocity_gradient_tensor(Particle current, std::vector<Particle>& neighbors, double dim, bool use_tensor_correction);
void compute_stress_rate(std::vector<Particle>& particles, double G, double dim, bool use_tensor_correction);

