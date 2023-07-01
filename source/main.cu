#include <iostream>
#include <string>
#include "../include/vec.cuh"
#include "../include/constants.h"
#include "../include/particle.cuh"
#include "../include/vtkwrite.h"
#include "../include/vtkgrid.h"



__global__ void compute_position(Particle* particles, size_t n, double dt){
    size_t i = blockDim.x * blockIdx.x + threadIdx.x;
    if(i < n){ //for every particle
        //Position x(t+dt)
        particles[i].pos = particles[i].pos + particles[i].vel * dt + particles[i].acc * dt *dt / 2;
        //velocity intermediate v(t+dt/2)
        particles[i].vel = particles[i].vel + particles[i].acc * dt / 2;

        //reflect bounary condition at y = -1
        if(particles[i].pos.z <= -1 || particles[i].pos.z >= 3)
            {particles[i].vel.z = -particles[i].vel.z;}
        
        if(particles[i].pos.x <= -1 || particles[i].pos.x >= 3)
           { particles[i].vel.x = -particles[i].vel.x;}
        
        if(particles[i].pos.y <= -1 || particles[i].pos.y >= 3)
           { particles[i].vel.y = -particles[i].vel.y;}
                
        // printf("Particle %d - Velocity: x = %f, y = %f, z = %f\n", i, particles[i].vel.x, particles[i].vel.y, particles[i].vel.z);

    }
}

//Heavyside function used in granular force computation
__host__ __device__ bool step(double x){
    return (x < 0) ? 0 : 1;
}

// -------------------New CODE----------------By Chirag----------------//
__global__ void compute_acceleration_sd(Particle* particles, size_t n, double dt, Vector3d g, Vector3d* force_d, Const c_h){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    double ke = 0, pe = 0;

    if(i < n){//acceleration a(t+dt)

        ke += 0.5 * dot(particles[i].vel, particles[i].vel) * c_h.M;

        Vector3d temp_force= Vector3d (0,0,0);

        for(size_t j = 0; j < n; j++){

            if (i==j){
                continue;
            }

            Vector3d x_ij = particles[i].pos - particles[j].pos;
            double x_norm = x_ij.norm();
            Vector3d x_unit = x_ij / x_norm;
            
            temp_force = c_h.K * x_unit * (c_h.sigma - x_norm) * step(c_h.sigma - x_norm) - 
                c_h.gamma * x_unit * dot(x_unit, (particles[i].vel - particles[j].vel)) * step(c_h.sigma - x_norm);//Granular forces(Spring-Dashpot) 
                
            printf("Particle %d - Temp Force: (%f, %f, %f)\n", i, temp_force.x, temp_force.y, temp_force.z);

            pe += 4.0 * c_h.eps * (pow(c_h.sigma / x_norm, 12) - pow(c_h.sigma / x_norm, 6));
        }

        particles[i].acc = temp_force/ c_h.M;


    }
}

// -------------------OLD CODE----------------By Praneeth----------------//

// __global__ void compute_acceleration_sd(Particle* particles, size_t n, double dt, Vector3d g, Vector3d* force_d, Const c_h){
//     int i = blockIdx.x * blockDim.x + threadIdx.x;
//     double ke = 0, pe = 0;
//     if(i < n){//acceleration a(t+dt)
//         ke += 0.5 * dot(particles[i].vel, particles[i].vel) * c_h.M;
//         for(size_t j = i + 1; j < n; j++){
//             Vector3d x_ij = particles[i].pos - particles[j].pos;
//             double x_norm = x_ij.norm();
//             Vector3d x_unit = x_ij / x_norm;
//             Vector3d f;
//             f = c_h.K * x_unit * (c_h.sigma - x_norm) * step(c_h.sigma - x_norm) 
//             - c_h.gamma * x_unit * dot(x_unit, (particles[i].vel - particles[j].vel)) * step(c_h.sigma - x_norm);//Granular forces(Spring-Dashpot) 
//             force_d[j] = force_d[j] - f;//reaction force
//             force_d[i] = force_d[i] + f;
//             pe += 4.0 * c_h.eps * (pow(c_h.sigma / x_norm, 12) - pow(c_h.sigma / x_norm, 6));
//         }
//         __syncthreads();
//         force_d[i] = force_d[i] + c_h.M * g;
//         particles[i].acc = force_d[i] / c_h.M;
//     }
//     // std::cout<<"Internal energy:"<<ke + pe<<std::endl;
// }

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
    size_t CR = ceil(cbrtf(n));
    double T = std::stod(argv[2]);
    double dt = std::stod(argv[3]);
    double gx = std::stod(argv[4]);
    double gy = std::stod(argv[5]);
    double gz = std::stod(argv[6]);

    Vector3d g(gx, gy, gz);

    //Create constants
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
        particles[i] = Particle(i + 1);
        // particles[i].initialize_grid(i);
    }

    //OR
    //read particles from files
    initialize_from_file(particles, "../initial_conditions.txt");

    // Create bounding box
    create_wireframe(-1.5, 3.5, -1.5, 3.5, -1.5, 3.5);

    //copy *particles to gpu
    size_t size = n * sizeof(Particle);
    cudaMalloc(&particles_device, size);
    cudaMemcpy(particles_device, particles, size, cudaMemcpyHostToDevice);

    int threadsPerBlock = 512;
    int blocksPerGrid = (n + threadsPerBlock - 1) / threadsPerBlock;

    int iter = 0;
    size_t size_v = n * sizeof(Vector3d);
    for(double t = 0; t<=T; t += dt){
        // std::cout<<"time:"<<t<<" "<<std::endl;
        
        Vector3d* force_d;
        cudaMalloc(&force_d, size_v);

        compute_position<<<blocksPerGrid, threadsPerBlock>>>(particles_device, n, dt);
        cudaDeviceSynchronize();

        compute_acceleration_sd<<<blocksPerGrid, threadsPerBlock>>>(particles_device, n, dt, g, force_d, c_h);
        cudaDeviceSynchronize();

        compute_velocity<<<blocksPerGrid, threadsPerBlock>>>(particles_device, n, dt);
        cudaDeviceSynchronize();
        
        cudaMemcpy(particles, particles_device, size, cudaMemcpyDeviceToHost);

        if(iter%10 == 0)
            
            write_vtk(particles, n, iter);
        iter++;

        cudaFree(force_d);
    }

    cudaFree(particles_device);
    delete[] particles;
}
