#pragma once
#include "Data_structs.hpp"
#include "Kernel.hpp"
#include "Hydro.hpp"
using namespace Data_structs;

void sph_leapfrog_step(
    std::vector<Particle>& particles,
    double dt,
    double dim,
    double k,
    double gamma,
    bool use_tensor_correction,
    bool use_wendland);

void density_update(
    std::vector<Particle>& particles,
    double dt,
    double dim,
    bool use_wendland,
    bool use_tensor_correction);