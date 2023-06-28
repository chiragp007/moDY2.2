#include "../include/particle.h"

Particle::Particle() {}

Particle::Particle(int id, std::mt19937& gen) : id(id),init_cell_index(-1) {

    std::normal_distribution<double> distribution(0.0, Const::std_dev*0.01);
    vel.set(distribution(gen), distribution(gen), distribution(gen));
    acc.set(0,0,0);
    pos = double(id - 1) * Const::imd * Vector3d(1, 0, 0);
    
}