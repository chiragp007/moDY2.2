#pragma once

#include <random>
#include "vec.h"
#include "constants.h"

class Particle{
    public:
        int id;
        Vector3d pos, vel, acc;
        int init_cell_index;
        int next_particle;

        Particle(); //Default Constructor
        Particle(int id, std::mt19937& gen); //Constructor assigns random velocities and id to the particles
};