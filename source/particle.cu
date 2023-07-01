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
    std::getline(file, line);
    std::getline(file, line);
    for(size_t i = 0; i < n; i++){
        file >> particles[i].pos.x >> particles[i].pos.y >> particles[i].pos.z;  
    }
    std::getline(file, line);
    std::getline(file, line);
    for(size_t i = 0; i < n; i++){
        file >> particles[i].vel.x >> particles[i].vel.y >> particles[i].vel.z;  
    }    

    std::cout<<"n "<<n<<std::endl;
    for(size_t i = 0; i < n; i++){
        std::cout<<particles[i].vel.x <<" "<<particles[i].vel.y << " "<< particles[i].vel.z<<std::endl; 
    } 

    file.close();
}


