#pragma once

#include "vec.cuh"
#include "constants.h"
#include <random>

class Particle{
    public:
        int id;
        Vector3d pos{}, vel{}, acc{};

        Particle() = default; //Default Constructor
        Particle(int id, std::mt19937& gen, Const c); //Constructor assigns random velocities and id to the particles
};