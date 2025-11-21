#include <iostream>
#include <vector>
#include <string>
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
int main(int argc, char* argv[]){
    std::string filename;
    int steps = 100;
    bool use_shepard = false;
    bool use_tensor_correction = false;
    bool use_consistent_shepard = false;
    bool use_wendland = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-f") {
            if (i + 1 < argc) {
                filename = argv[++i];
            } else {
                std::cerr << "Error: -f option requires a filename." << std::endl;
                return 1;
            }
        } else if (arg == "-N") {
            if (i + 1 < argc) {
                try {
                    steps = std::stoi(argv[++i]);
                } catch (const std::invalid_argument& e) {
                    std::cerr << "Error: Invalid number for -N option." << std::endl;
                    return 1;
                }
            } else {
                std::cerr << "Error: -N option requires a number." << std::endl;
                return 1;
            }
        } else if (arg == "-S") {
            use_shepard = true;
        } else if (arg == "-T") {
            use_tensor_correction = true;
        } else if (arg == "-CS") {
            // Consistent Shepard interpolation flag
            // This could be handled as needed in the simulation loop
            // For now, we just print a message
            std::cout << "Consistent Shepard interpolation will be applied." << std::endl;
            use_consistent_shepard = true;
        } else if (arg == "-W") {
            use_wendland = true;
        }
         else if (arg == "-h" || arg == "--help") {
            std::cout << "Usage: " << argv[0] << " -f <filename> [-N <steps>] [-S] [-T]" << std::endl;
            std::cout << "  -f <filename>    Input HDF5 file with initial particle data" << std::endl;
            std::cout << "  -N <steps>       Number of simulation steps (default: 100)" << std::endl;
            std::cout << "  -S               Use Shepard correction" << std::endl;
            std::cout << "  -T               Use Tensor correction" << std::endl;
            std::cout << "  -CS              Use Consistent Shepard interpolation" << std::endl;
            std::cout << "  -W               Use Wendland kernel for density computation" << std::endl;
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
    
    // Simulation parameters
    double CFL = 0.3;
    double k = 1.0;
    double gamma = 1.4;
    compute_density(particles, 3, false, use_wendland);
    if(use_shepard){
        shepard_correction(particles, 3, use_wendland);
        compute_density(particles, 3, true, use_wendland);
    }
    if(use_consistent_shepard){
        consistent_shepard_interpolation(particles, 3, use_wendland);
        compute_density(particles, 3, true, use_wendland);
    }
    if(use_tensor_correction){
        tensor_correction(particles, 3, use_wendland);
    }
    compute_pressure(particles, k, gamma);
    compute_sound_speed(particles, gamma);
    compute_acceleration(particles, 3, use_tensor_correction, use_wendland);
    double dt = CFL * particles[0].smoothing_length / max_sound_speed(particles);
    // Main simulation loop
    for(int i = 0; i < steps; ++i) {
        double max_cs = max_sound_speed(particles);
        sph_leapfrog_step(particles, dt, 3, k, gamma, use_tensor_correction, use_wendland);
        std::cout << "Completed step " << i+1 << "/" << steps << std::endl;
        // Update time step based on CFL condition
        if (max_cs > 0) {
            dt = CFL * particles[0].smoothing_length / max_cs; // CFL condition update
            std::cout << "Updated time step to " << dt << std::endl;
        }
    }
    write_output_particles("output_particles.csv", particles);
    std::cout << "SPH Simulation Ended" << std::endl;
    return 0;
}