#include "../include/particle.cuh"

Particle::Particle(int id, std::mt19937& gen, Const c) : id(id) {
    std::normal_distribution<double> distribution(0.0, c.std_dev);
    vel.set(distribution(gen), distribution(gen), distribution(gen));
    acc.set(0,0,0);
    pos.set((double(id - 1) * c.imd * 1), 0, 0);
}