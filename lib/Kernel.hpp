#pragma once
#include "Data_structs.cpp"
#include <vector>
#include <cmath>
#include <stdexcept>
using namespace Data_structs;

double W(double r, double h, int dim);
double derivW(double r, double h, int dim);
Vector gradW(Vector pos_i, Vector pos_j, double h, int dim);
void consistent_shepard_interpolation(std::vector<Particle>& particles, double h, int dim);
void tensor_Correction(std::vector<Particle>& particles, double h, int dim);