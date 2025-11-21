#pragma once
#include "Data_structs.hpp"

using namespace Data_structs;

double W(double r, double h, int dim);
double wendland_W(double r, double h, int dim);
double wendland_W_deriv(double r, double h, int dim);
double derivW(double r, double h, int dim);
Vector gradW(Vector pos_i, Vector pos_j, double h, int dim, bool use_wendland);
void shepard_correction(std::vector<Particle>& particles, int dim, bool use_wendland);
void tensor_correction(std::vector<Particle>& particles, int dim, bool use_wendland);
void consistent_shepard_interpolation(std::vector<Particle>& particles, int dim, bool use_wendland);