
#include <fstream>
#include <iostream>
#include <string>
#include "particle.h"

void write_vtk(Particle* particles, size_t n, int iter, size_t grid_size_x, size_t grid_size_y, size_t grid_size_z) {
    std::string s = "../output/t" + std::to_string(iter) + ".vtk";
    std::ofstream file(s);

    // Write the header information
    file << "# vtk DataFile Version 2.0" << std::endl;
    file << "Particle Simulation Data" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;

    // Write the points information
    file << "POINTS " << n << " double" << std::endl;
    for (size_t i = 0; i < n; ++i) {
        file << particles[i].pos.x << " " << particles[i].pos.y << " " << particles[i].pos.z << std::endl;
    }

    // Write the grid information
    file << "CELLS " << n << " " << 2 * n << std::endl;
    for (size_t i = 0; i < n; ++i) {
        file << "1 " << i << std::endl;
    }

    // Write the cell types
    file << "CELL_TYPES " << n << std::endl;
    for (size_t i = 0; i < n; ++i) {
        file << "1" << std::endl;
    }

    // Write the point data (velocity)
    file << "POINT_DATA " << n << std::endl;
    file << "VECTORS velocity float" << std::endl;
    for (size_t i = 0; i < n; ++i) {
        file << particles[i].vel.x << " " << particles[i].vel.y << " " << particles[i].vel.z << std::endl;
    }

    // Write the grid dimensions as metadata
    file << std::endl;
    file << "METADATA" << std::endl;
    file << "grid_size_x " << grid_size_x << std::endl;
    file << "grid_size_y " << grid_size_y << std::endl;
    file << "grid_size_z " << grid_size_z << std::endl;

    // Close the file
    file.close();
}
