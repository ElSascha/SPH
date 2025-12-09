#pragma once
#include "Data_structs.hpp"
#include <vector>
#include <unordered_map>

using namespace Data_structs;

// Spatial hash key from 3D cell indices
inline size_t hash_cell(int ix, int iy, int iz) {
    // Large primes for hashing
    return static_cast<size_t>(ix * 73856093) ^ 
           static_cast<size_t>(iy * 19349663) ^ 
           static_cast<size_t>(iz * 83492791);
}

class NeighborSearch {
public:
    double cell_size;
    std::unordered_map<size_t, std::vector<size_t>> cells;
    
    NeighborSearch(double h) : cell_size(2.0 * h) {}
    
    void build(const std::vector<Particle>& particles) {
        cells.clear();
        for (size_t i = 0; i < particles.size(); ++i) {
            int ix = static_cast<int>(std::floor(particles[i].position.x / cell_size));
            int iy = static_cast<int>(std::floor(particles[i].position.y / cell_size));
            int iz = static_cast<int>(std::floor(particles[i].position.z / cell_size));
            cells[hash_cell(ix, iy, iz)].push_back(i);
        }
    }
    
    // Get neighbor indices for particle at position pos within radius 2h
    template<typename Callback>
    void for_each_neighbor(const Vector& pos, double h, Callback&& callback) const {
        double search_radius = 2.0 * h;
        int ix0 = static_cast<int>(std::floor((pos.x - search_radius) / cell_size));
        int iy0 = static_cast<int>(std::floor((pos.y - search_radius) / cell_size));
        int iz0 = static_cast<int>(std::floor((pos.z - search_radius) / cell_size));
        int ix1 = static_cast<int>(std::floor((pos.x + search_radius) / cell_size));
        int iy1 = static_cast<int>(std::floor((pos.y + search_radius) / cell_size));
        int iz1 = static_cast<int>(std::floor((pos.z + search_radius) / cell_size));
        
        for (int ix = ix0; ix <= ix1; ++ix) {
            for (int iy = iy0; iy <= iy1; ++iy) {
                for (int iz = iz0; iz <= iz1; ++iz) {
                    auto it = cells.find(hash_cell(ix, iy, iz));
                    if (it != cells.end()) {
                        for (size_t j : it->second) {
                            callback(j);
                        }
                    }
                }
            }
        }
    }
};
