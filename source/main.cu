#include <iostream>
#include <string>
#include <random>
#include <cuda_runtime_api.h>
#include <cuda_runtime.h>

#include "../include/vec.cuh"
#include "../include/constants.h"
#include "../include/particle.cuh"
#include "../include/vtkwrite.h"



__device__ double energy_d = 0;

//Heavyside function used in granular force computation
__host__ __device__ bool step(double x){
    return (x < 0) ? 0 : 1;
}

__global__ void compute_position(Particle* particles, size_t n, double dt){
    size_t i = blockDim.x * blockIdx.x + threadIdx.x;
    if(i < n){ //for every particle
        //Position x(t+dt)
        particles[i].pos = particles[i].pos + particles[i].vel * dt + particles[i].acc * dt *dt / 2;

        //velocity intermediate v(t+dt/2)
        particles[i].vel = particles[i].vel + particles[i].acc * dt / 2;
    }
}



__global__ void compute_acceleration(Particle* particles, size_t n, double dt, Const c_d,Vector3d* force_d){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    double ke = 0.0f, pe = 0.0f;
    
    if(i < n){     
        ke += 0.5 * dot(particles[i].vel, particles[i].vel) * c_d.M;
        force_d[i] = Vector3d(0, 0, 0);

        for(int j =0; j < n; j++){

            if (i==j){
                continue;
            }

            Vector3d x_ij = particles[i].pos - particles[j].pos;
            double x_norm = x_ij.norm();

            Vector3d temp1= (x_ij / pow(x_norm, 2));
            
            force_d[i] = force_d[i] + 24.0 * c_d.eps * pow(c_d.sigma / x_norm, 6) * (2 * pow(c_d.sigma / x_norm, 6) - 1)*temp1 ;//Atomic forces(Lennard-Jones)
           // force_d[i] = force_d[i] + (c_d.K * x_ij / x_norm) * (c_d.sigma - x_norm) * step(c_d.sigma - x_norm);//Granular forces(Spring-Dashpot)
            //force_d[j] = force_d[j] - force_d[i];//reaction force       

            pe += 4.0 * c_d.eps * (pow(c_d.sigma / x_norm, 12) - pow(c_d.sigma / x_norm, 6));   
            
        }
        
        energy_d = ke + pe;

        particles[i].acc = force_d[i] / c_d.M;
                       
    }
   
}


__global__ void compute_velocity(Particle* particles, size_t n, double dt){
    size_t i = blockDim.x * blockIdx.x + threadIdx.x;
    if(i < n){
        //velocity v(t+dt)
        particles[i].vel = particles[i].vel + particles[i].acc * dt / 2;
    }
}

int main(int argc, char* argv[]){
    //command line arguments
    size_t n = std::stoi(argv[1]);
    double T = std::stod(argv[2]);
    double dt = std::stod(argv[3]);
    
    const Const c_h("../constants.txt");
    std::cout<<c_h.M<<" "<<c_h.R<<std::endl;
    //Initialize random seed
    std::random_device rd;
    std::mt19937 gen(rd());
    
    //Dynamically allocate n particles
    Particle* particles = new Particle[n];
    Particle* particles_device; 
    for(size_t i = 0; i < n; i++){
        // Initialize the particles
        particles[i] = Particle(i + 1, gen, c_h);
    }

    size_t size = n * sizeof(Particle);
    cudaMalloc(&particles_device, size);
    cudaMemcpy(particles_device, particles, size, cudaMemcpyHostToDevice);

    int threadsPerBlock = 512;
    int blocksPerGrid = (n + threadsPerBlock - 1) / threadsPerBlock;

    int iter = 0;
    double energy = 0;

    Vector3d force_h;
    Vector3d* force_d;
    cudaMalloc(&force_d, sizeof(Vector3d));
    
    for(double t = 0; t<=T; t += dt){

        std::cout<<"Time:"<<t<<" ";
        // Compute force after every time step using leapfrog integration scheme
       
        compute_position<<<blocksPerGrid, threadsPerBlock>>>(particles_device, n, dt);
        compute_acceleration<<<blocksPerGrid, threadsPerBlock>>>(particles_device, n, dt, c_h,force_d);

        compute_velocity<<<blocksPerGrid, threadsPerBlock>>>(particles_device, n, dt);
        
        
        // this is needed since we write it after every dt
        cudaMemcpy(particles, particles_device, size, cudaMemcpyDeviceToHost);

        if(iter%100 == 0)
            write_vtk(particles, n, iter);

        cudaMemcpyFromSymbol(&energy, energy_d, sizeof(double));
        
        printf("Internal energy:%f\n",energy);
        
        cudaFree(force_d);
        iter++;
    }

    cudaFree(particles_device);
    delete[] particles;
    //delete[] force_d;
//    cudaFree(force_d)
}