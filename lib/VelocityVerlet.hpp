#pragma once
#include <vector>
#include "Data_structs.cpp"
#include "Kernel.hpp"
#include "Hydro.hpp"
using namespace Data_structs;

void velocity_verlet_step(std::vector<Particle>& particles, double dt, double h, double dim, bool use_shepard, double k, double gamma);