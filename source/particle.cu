#include "../include/particle.cuh"
#include <fstream>
#include <iostream>

Particle::Particle() {}

Particle::Particle(int id) : id(id) {
}

void Particle::initialize_grid(int i, size_t CR, Const c){
    this->pos.x = ((i / CR) % CR) * c.imd;
    this->pos.y = (i / (CR * CR)) * c.imd;
    this->pos.z = (i % CR) * c.imd;   
}

void Particle::initialize_random(){

}

void initialize_from_file(Particle* particles, std::string filepath){
    std::ifstream file;
    file.open(filepath);

    size_t n;
    std::string line;
    
    file >> n;

    //read radius
    std::getline(file, line);
    std::getline(file, line);
    for(size_t i = 0; i < n; i++){
        file >> particles[i].radius;  
        // M is 5
        particles[i].moi = 5* particles[i].radius * particles[i].radius;
    }

    //Pos
    std::getline(file, line);
    std::getline(file, line);
    for(size_t i = 0; i < n; i++){
        file >> particles[i].pos.x >> particles[i].pos.y >> particles[i].pos.z;  
    }

    //Vel
    std::getline(file, line);
    std::getline(file, line);
    for(size_t i = 0; i < n; i++){
        file >> particles[i].vel.x >> particles[i].vel.y >> particles[i].vel.z;  
    }    

    // Init_rotation
    std::getline(file, line);
    std::getline(file, line);
    for(size_t i = 0; i < n; i++){
        file >> particles[i].ini_rotate.x >> particles[i].ini_rotate.y >> particles[i].ini_rotate.z;  
    }


    // Ang_vel
    std::getline(file, line);
    std::getline(file, line);
    for(size_t i = 0; i < n; i++){
        file >> particles[i].ang_vel.x >> particles[i].ang_vel.y >> particles[i].ang_vel.z;  
    }

    file.close();
}


