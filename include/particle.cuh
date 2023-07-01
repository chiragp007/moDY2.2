#pragma once

#include <random>
#include <cmath>
#include "vec.cuh"
#include "constants.h"

extern size_t N;
extern size_t CR; 

class Particle{
    public:
        int id;
        Vector3d pos, vel, acc; 

        Particle(); //Default Constructor
        Particle(int id); //Constructor assigns random velocities and id to the particles
        
        void initialize_grid(int i, size_t CR, Const c);
        void initialize_random();
};

void initialize_from_file(Particle* particles, std::string filename);