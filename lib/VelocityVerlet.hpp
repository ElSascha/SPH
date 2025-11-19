#pragma once
#include "Data_structs.hpp"
#include "Kernel.hpp"
#include "Hydro.hpp"
using namespace Data_structs;

void velocity_verlet_step(std::vector<Particle>& particles, double dt, double dim, bool use_shepard,bool use_tensor_correction, double k, double gamma, bool use_wendland);