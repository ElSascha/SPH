#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <limits>
#include <cmath>
#include <H5Cpp.h>
#include "../lib/Data_structs.hpp"
#include "../lib/VelocityVerlet.hpp"
#include "../lib/Hydro.hpp"

using namespace Data_structs;

std::vector<Particle> load_initial_particles(const std::string& filename) {
    std::vector<Particle> particles;

    try {
        H5::H5File file(filename, H5F_ACC_RDONLY);

        // --- Positionen ---
        H5::DataSet ds_x = file.openDataSet("positions");
        H5::DataSpace space_x = ds_x.getSpace();
        hsize_t dims_x[2];
        space_x.getSimpleExtentDims(dims_x);
        size_t N = dims_x[0];
        std::vector<double> pos_data(N * 3);
        ds_x.read(pos_data.data(), H5::PredType::NATIVE_DOUBLE);

        // --- Geschwindigkeiten ---
        H5::DataSet ds_v = file.openDataSet("velocities");
        std::vector<double> vel_data(N * 3);
        ds_v.read(vel_data.data(), H5::PredType::NATIVE_DOUBLE);

        // --- Massen ---
        H5::DataSet ds_m = file.openDataSet("mass");
        std::vector<double> mass_data(N);
        ds_m.read(mass_data.data(), H5::PredType::NATIVE_DOUBLE);

        // --- Smoothing Lengths ---
        H5::DataSet ds_h = file.openDataSet("smoothing_length");
        std::vector<double> sml_data(N);
        ds_h.read(sml_data.data(), H5::PredType::NATIVE_DOUBLE);

        // --- Partikel zusammensetzen ---
        particles.reserve(N);
        for (size_t i = 0; i < N; ++i) {
            Particle p;
            p.position = Vector(pos_data[i*3 + 0],
                                pos_data[i*3 + 1],
                                pos_data[i*3 + 2]);
            p.velocity = Vector(vel_data[i*3 + 0],
                                vel_data[i*3 + 1],
                                vel_data[i*3 + 2]);
            p.mass = mass_data[i];
            p.smoothing_length = sml_data[i];

            // Initialize default physical values
            p.acceleration = Vector(0, 0, 0);
            p.density = 0.0;
            p.rho_pred = 0.0;
            p.drho_dt = 0.0;
            p.pressure = 0.0;
            p.sound_speed = 0.0;
            p.shepard = 0.0;

            particles.push_back(p);
        }

        std::cout << "Loaded " << N << " particles from " << filename << std::endl;
    }
    catch (const H5::FileIException& e) {
        std::cerr << "HDF5 File Error: " << e.getDetailMsg() << std::endl;
    }
    catch (const H5::DataSetIException& e) {
        std::cerr << "HDF5 Dataset Error: " << e.getDetailMsg() << std::endl;
    }
    catch (const H5::DataSpaceIException& e) {
        std::cerr << "HDF5 Dataspace Error: " << e.getDetailMsg() << std::endl;
    }

    return particles;
}

void write_output_particles(const char* filename, const std::vector<Particle>& particles){
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }
    // Write header
    file << "x,y,z,vx,vy,vz,mass,smoothing_length,density,pressure,sound_speed\n";
    for (const auto& p : particles) {
        file << p.position.x << "," << p.position.y << "," << p.position.z << ","
             << p.velocity.x << "," << p.velocity.y << "," << p.velocity.z << ","
             << p.mass << "," << p.smoothing_length << "," << p.density << "," << p.pressure << "," << p.sound_speed << "\n";
    }
    file.close();
}

double max_sound_speed(const std::vector<Particle>& particles){
    double max_cs = 0.0;
    for(const auto& p : particles){
        if(p.sound_speed > max_cs){
            max_cs = p.sound_speed;
        }
    }
    return max_cs;
}

// Longitudinal wave speed for solids: c_L = sqrt(c_bulk^2 + 4/3 * G / rho)
double max_longitudinal_speed(const std::vector<Particle>& particles, double G){
    double max_c = 0.0;
    for(const auto& p : particles){
        if(p.density <= 0.0) continue;
        double c_bulk = p.sound_speed;
        double cL_sq = c_bulk * c_bulk + (4.0/3.0) * G / p.density;
        if(cL_sq > 0.0){
            double cL = std::sqrt(cL_sq);
            if(cL > max_c) max_c = cL;
        }
    }
    return max_c;
}

double min_smoothing_length(const std::vector<Particle>& particles){
    double h_min = std::numeric_limits<double>::infinity();
    for(const auto& p : particles){
        if(p.smoothing_length > 0.0 && p.smoothing_length < h_min){
            h_min = p.smoothing_length;
        }
    }
    return h_min;
}
int main(int argc, char* argv[]){
    std::string filename;
    int total_time = 1;
    double time_step = 0.01;
    bool use_shepard = false;
    bool use_tensor_correction = false;
    bool use_consistent_shepard = false;
    bool test = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-f") {
            if (i + 1 < argc) {
                filename = argv[++i];
            } else {
                std::cerr << "Error: -f option requires a filename." << std::endl;
                return 1;
            }
        } else if (arg == "-total_time") {
            if (i + 1 < argc) {
                try {
                    total_time = std::stoi(argv[++i]);
                } catch (const std::invalid_argument& e) {
                    std::cerr << "Error: Invalid number for -N option." << std::endl;
                    return 1;
                }
            } else {
                std::cerr << "Error: -N option requires a number." << std::endl;
                return 1;
            }
        }else if(arg == "-time_step"){
            if (i + 1 < argc) {
                try {
                    time_step = std::stod(argv[++i]);
                } catch (const std::invalid_argument& e) {
                    std::cerr << "Error: Invalid number for -time_step option." << std::endl;
                    return 1;
                }
            } else {
                std::cerr << "Error: -time_step option requires a number." << std::endl;
                return 1;
            }
        }else if (arg == "-S") {
            use_shepard = true;
            std::cout << "Shepard correction enabled." << std::endl;
        } else if (arg == "-T") {
            use_tensor_correction = true;
            std::cout << "Tensor correction enabled." << std::endl;
        } else if (arg == "-CS") {
            // Consistent Shepard interpolation flag
            // This could be handled as needed in the simulation loop
            // For now, we just print a message
            std::cout << "Consistent Shepard interpolation will be applied." << std::endl;
            use_consistent_shepard = true;
        } else if (arg == "-test"){
            test = true;
            std::cout << "Test mode enabled." << std::endl;
        }
         else if (arg == "-h" || arg == "--help") {
            std::cout << "Usage: " << argv[0] << " -f <filename> [-N <steps>] [-S] [-T]" << std::endl;
            std::cout << "  -f <filename>    Input HDF5 file with initial particle data" << std::endl;;
            std::cout << "  -S               Use Shepard correction" << std::endl;
            std::cout << "  -T               Use Tensor correction" << std::endl;
            std::cout << "  -CS              Use Consistent Shepard interpolation" << std::endl;
            std::cout << "  -time_step <dt>  Time step size (default: 0.01)" << std::endl;
            std::cout << "  -test            Run in test mode (output initial state and exit)" << std::endl;
            std::cout << "  -h, --help       Show this help message" << std::endl;

            return 0;
        }
        else {
            std::cerr << "Error: Unknown option " << arg << std::endl;
            std::cerr << "Usage: " << argv[0] << " -f <filename> [-N <steps>] [-S] [-T]" << std::endl;
            return 1;
        }
    }

    if (filename.empty()) {
        std::cerr << "Error: No input file specified." << std::endl;
        std::cerr << "Usage: " << argv[0] << " -f <filename> [-N <steps>] [-S] [-T]" << std::endl;
        return 1;
    }

    std::cout << "SPH Simulation Started" << std::endl;
    // Initialize particles
    std::vector<Particle> particles = load_initial_particles(filename.c_str());
    if (particles.empty()) {
        std::cerr << "Error: Failed to load particles. Exiting." << std::endl;
        return 1;
    }
    
    // Simulation parameters for Basalt Murnaghan EOS
    double CFL = 0.3;               
    double K_0 = 40e9; // Bulk modulus in Pascals
    double K_0_deriv = 4.0; // Derivative of bulk modulus
    double E = 50e9; // Youngs modulus in Pascals
    double nu = 0.25; // Poisson's ratio
    double G = E / (2.0 * (1.0 + nu)); // Shear modulus
    
    compute_density(particles, 3, false);
    if(use_shepard){
        shepard_correction(particles, 3);
        compute_density(particles, 3, true);
    }
    if(use_consistent_shepard){
        consistent_shepard_interpolation(particles, 3);
        compute_density(particles, 3, true);
    }
    if(use_tensor_correction){
        tensor_correction(particles, 3);
    }
    // Initialize rho_0 to current density to ensure zero pressure at start
    for(auto& p : particles){
        p.rho_0 = p.density;
    }

    compute_pressure(particles, K_0, K_0_deriv);
    compute_sound_speed(particles, K_0, K_0_deriv);
    compute_stress_rate(particles, G, 3, use_tensor_correction);
    // Artificial viscosity parameters for basalt
    double alpha_visc = 1.0;   // Linear viscosity (shear)
    double beta_visc = 2.0;    // Quadratic viscosity (shock)
    double epsilon_visc = 0.01; // Singularity prevention
    compute_acceleration(particles, 3, use_tensor_correction, alpha_visc, beta_visc, epsilon_visc);
    if(test){
        write_output_particles("output_particles.csv", particles);
        std::cout << "Test output written. Exiting." << std::endl;
        return 0;
    }
    double dt = 0;
    
    // Main simulation loop
    double current_time = 0.0;
    double next_output_time = 0.0;
    while(total_time > current_time){
        double max_cL = max_longitudinal_speed(particles, G);
        double h_min = min_smoothing_length(particles);
        
        // Find maximum velocity for viscosity CFL
        double max_v = 0.0;
        for(const auto& p : particles){
            double v = p.velocity.norm();
            if(v > max_v) max_v = v;
        }

        if (max_cL > 0.0 && std::isfinite(h_min)) {
            // CFL condition including artificial viscosity:
            // dt <= CFL * h / (c_L + alpha*c + beta*h*|div v|_max)
            // Simplified: use max velocity as proxy for h*|div v|
            double signal_speed = max_cL + alpha_visc * max_cL + beta_visc * max_v;
            dt = CFL * h_min / signal_speed;
            std::cout << "Updated time step to " << dt << " (h_min=" << h_min << ", c_L=" << max_cL << ")" << std::endl;
        }
        else {
            std::cerr << "Warning: Wave speed or smoothing length invalid. Time step not updated." << std::endl;
        }
        sph_leapfrog_step(
            particles,
            dt,
            3,
            K_0,
            K_0_deriv,
            G,
            use_tensor_correction,
            use_shepard,
            use_consistent_shepard);
        std::cout << "Completed time: " << current_time + dt << " / " << total_time << std::endl;
        // Update time step based on CFL condition
        // write output when time_step is reached
        current_time += dt;
        if(current_time >= next_output_time) {
            write_output_particles(("data_test/output_particles_step_" + std::to_string(static_cast<int>(next_output_time / time_step)) + ".csv").c_str(), particles);
            next_output_time += time_step;
        }

    }
    std::cout << "SPH Simulation Ended" << std::endl;
    return 0;
}