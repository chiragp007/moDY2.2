#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include "particle.cuh"

void write_vtk(Particle* particles, size_t n, int iter) {

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

    // Write the point data (velocity)
    file<<"CELLS 0 0"<<std::endl;
    file<<"CELL_TYPES 0"<<std::endl;
    file << "POINT_DATA " << n << std::endl;
    file << "VECTORS velocity float" << std::endl;
    for (size_t i = 0; i < n; ++i) {
    file << particles[i].vel.x << " " << particles[i].vel.y << " " << particles[i].vel.z << std::endl;
    }

    // Close the file
    file.close();
}