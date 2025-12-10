#pragma once
#include "Data_structs.hpp"
#include "Kernel.hpp"
#include "Hydro.hpp"
using namespace Data_structs;

void sph_leapfrog_step(
    std::vector<Particle>& particles,
    double dt,
    double dim,
    double K_0,
    double K_0_deriv,
    double G,
    double alpha_visc,
    double beta_visc,
    double epsilon_visc,
    bool use_tensor_correction,
    bool use_shepard,
    bool use_consistent_shepard,
    bool integrate_density
    );
