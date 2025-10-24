#pragma once
#include "Data_structs.hpp"
#include <vector>
#include <cmath>
#include <stdexcept>
using namespace Data_structs;

double W(double r, double h, int dim);
double derivW(double r, double h, int dim);
Vector gradW(Vector pos_i, Vector pos_j, double h, int dim);
void shepard_correction(std::vector<Particle>& particles, double h, int dim);
void tensor_correction(std::vector<Particle>& particles, double h, int dim);