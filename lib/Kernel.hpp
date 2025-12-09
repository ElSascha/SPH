#pragma once
#include "Data_structs.hpp"

using namespace Data_structs;

double W(double r, double h, int dim);
double W_deriv(double r, double h, int dim);
Vector gradW(Vector pos_i, Vector pos_j, double h, int dim);

// Combined kernel + gradient computation (faster when you need both)
void W_and_gradW(double r, double h, const Vector& r_vec, double& W_out, Vector& gradW_out);

void shepard_correction(std::vector<Particle>& particles, int dim);
void tensor_correction(std::vector<Particle>& particles, int dim);
void consistent_shepard_interpolation(std::vector<Particle>& particles, int dim);
