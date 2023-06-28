#include <iostream>
#include <string>
#include "../include/vec.h"
#include "../include/constants.h"
#include "../include/particle.h"
#include "../include/vtkwrite.h"


struct extended_list {
    int integer;
    int ID;

};


void compute_position(Particle* particles, size_t n, double dt, int  grid_size_x, int grid_size_y, int grid_size_z, extended_list* list_cells, extended_list* list_particles,int box_width, int no_cells,double spacing) {


    for (size_t i = 0; i < n; i++) { // For every particle
        // Position x(t+dt)
        particles[i].pos = particles[i].pos + particles[i].vel * dt + particles[i].acc * dt *dt / 2;

        // Velocity intermediate v(t+dt/2)
        particles[i].vel = particles[i].vel + particles[i].acc * dt / 2;

        if (particles[i].pos.x < 0)
            particles[i].pos.x += box_width;
        if (particles[i].pos.y < 0)
            particles[i].pos.y += box_width;
        if (particles[i].pos.z < 0)
            particles[i].pos.z += box_width;

        if (particles[i].pos.x >= box_width)
            particles[i].pos.x -= box_width;
        if (particles[i].pos.y >= box_width)
            particles[i].pos.y -= box_width;
        if (particles[i].pos.z >= box_width)
            particles[i].pos.z -= box_width;

        // while (particles[i].pos.x < 0)
        //     particles[i].pos.x += box_width;
        // while (particles[i].pos.y < 0)
        //     particles[i].pos.y += box_width;
        // while (particles[i].pos.z < 0)
        //     particles[i].pos.z += box_width;

        // while (particles[i].pos.x >= box_width)
        //     particles[i].pos.x -= box_width;
        // while (particles[i].pos.y >= box_width)
        //     particles[i].pos.y -= box_width;
        // while (particles[i].pos.z >= box_width)
        //     particles[i].pos.z -= box_width;

        // std::cout<<particles[i].pos.x<<" "<<particles[i].pos.y<<" "<<particles[i].pos.z<<std::endl;

        // Find the cell indices for the particle
        int cellX = particles[i].pos.x / spacing ;
        int cellY = particles[i].pos.y / spacing;
        int cellZ = particles[i].pos.z / spacing ;

        // int cell_index= cellX*grid_size_x*grid_size_x+ cellY*grid_size_x+ cellZ;
        // int cell_index = cellX + grid_size_x * (cellY + grid_size_y * cellZ);

        int cell_index = cellX + cellY * no_cells + cellZ * no_cells* no_cells;

        // std::cout << "Particle " << i << " - Position: (" << particles[i].pos.x << ", " << particles[i].pos.y << ", " << particles[i].pos.z << ") - Cell Index: " << cell_index << std::endl;

        if (cell_index < 0 || cell_index > no_cells*no_cells*no_cells) {
            std::cout << " Wrong Index"<<std::endl;
        } 

        // particles[i].next_particle = list_cells[cell_index].integer;

            // list_particles[i].integer = list_cells[cell_index].integer;

        // list_cells[cell_index].integer = list_particles[i].ID;
        // if (cell_index != particles[i].init_cell_index ){

        //     // std::cout << particles[i].init_cell_index <<" IF   "<<std::endl;

        //     particles[i].next_particle = list_cells[cell_index].integer;

        //     // list_particles[i].integer = list_cells[cell_index].integer;

        //     list_cells[cell_index].integer = list_particles[i].ID;
            
        //     particles[i].init_cell_index=cell_index; 

        // }
    }
}


// void compute_acceleration(Particle* particles, size_t n, double dt,int  grid_size_x, int grid_size_y, int grid_size_z, int spacing,extended_list* list_cells, extended_list* list_particles){
//     Vector3d* force = new Vector3d[n];
    
//     double ke = 0, pe = 0;
//     int no_cells=20;

//     for(size_t i = 0; i < n; i++){ //acceleration a(t+dt)
//         ke += 0.5 * dot(particles[i].vel, particles[i].vel) * Const::M;
        
//         int cellX = particles[i].pos.x / spacing;
//         int cellY = particles[i].pos.y / spacing;
//         int cellZ = particles[i].pos.z / spacing;

//         int cell_index = cellX + cellY * no_cells + cellZ * no_cells* no_cells;

//         int particleIndex = list_cells[cell_index].integer;
        
//          while (particleIndex != -1) {
//             // Calculate the force between particles in the same cell
//             if (particleIndex != i) {
//                 std::cout<<"Force "<<particleIndex<<std::endl;
//                 Vector3d x_ij = particles[i].pos - particles[particleIndex].pos;
                
//                 // if (x_ij.x>grid_size_x/2){
//                 //     x_ij.x -= grid_size_x;
//                 // }

//                 // if (x_ij.y>grid_size_y/2){
//                 //     x_ij.y -= grid_size_y;
//                 // }

//                 // if (x_ij.z>grid_size_z/2){
//                 //     x_ij.z -= grid_size_z;
//                 // }


//                 double x_norm = x_ij.norm();
//                 if (x_norm < 1.5) {
//                     force[i] = force[i] + 24.0 * Const::eps * pow(Const::sigma / x_norm, 6) * (2 * pow(Const::sigma / x_norm, 6) - 1) * (x_ij / pow(x_norm, 2));
//                     pe += 4.0 * Const::eps * (pow(Const::sigma / x_norm, 12) - pow(Const::sigma / x_norm, 6));
//                 }
//             }
//             particleIndex = list_particles[particleIndex].integer;
//         }

//         particles[i].acc = force[i] / Const::M;
//     }

//     std::cout << "Internal energy:" << ke + pe << std::endl;
//     delete[] force;
// }


void compute_acceleration(Particle* particles, size_t n, double dt,int  grid_size_x, int grid_size_y, int grid_size_z, int spacing,extended_list* list_cells, extended_list* list_particles){
    Vector3d* force = new Vector3d[n];
    int box_width=20;

    double ke = 0, pe = 0;
    int no_cells=20;
    for(size_t i = 0; i < n; i++){ //acceleration a(t+dt)
        ke += 0.5 * dot(particles[i].vel, particles[i].vel) * Const::M;
        
        int cellX = particles[i].pos.x / spacing;
        int cellY = particles[i].pos.y / spacing;
        int cellZ = particles[i].pos.z / spacing;

        int cell_index = cellX + cellY * no_cells + cellZ * no_cells* no_cells;

        int particleIndex = list_cells[cell_index].integer;
    
        for(size_t j = i + 1; j < n; j++){

            Vector3d x_ij = particles[i].pos - particles[j].pos;

            // std::cout << "BEFORE x_ij: (" << x_ij.x << ", " << x_ij.y << ", " << x_ij.z << ")" << std::endl;  

            // if (x_ij.x < 0)
            //     {x_ij.x += box_width;}
            
            // if (x_ij.x >= box_width)
            //     {x_ij.x -= box_width;}

            // if (x_ij.y < 0)
            //     {x_ij.y += box_width;}
            
            // if (x_ij.y >= box_width)
            //     {x_ij.y -= box_width;}

            // if (x_ij.z < 0)
            //     {x_ij.z += box_width;}
            
            // if (x_ij.z >= box_width)
            //     {x_ij.z -= box_width;}


            // std::cout << "AFTER x_ij: (" << x_ij.x << ", " << x_ij.y << ", " << x_ij.z << ")" << std::endl;      

            double x_norm = x_ij.norm();
            if (x_norm<2.5){
             force[i] = force[i] + 24.0 * Const::eps * pow(Const::sigma / x_norm, 6) * (2 * pow(Const::sigma / x_norm, 6) - 1) * (x_ij / pow(x_norm, 2));//Atomic forces(Lennard-Jones)
             pe += 4.0 * Const::eps * (pow(Const::sigma / x_norm, 12) - pow(Const::sigma / x_norm, 6));

            }
        }
        particles[i].acc = force[i] / Const::M;
    }
    std::cout<<"Internal energy:"<<ke + pe<<std::endl;
    delete[] force;
}


void compute_velocity(Particle* particles, size_t n, double dt){
    for(size_t i = 0; i < n; i++){
        //velocity v(t+dt)
        particles[i].vel = particles[i].vel + particles[i].acc * dt / 2;
    }
}

double Const::M, Const::R, Const::K, Const::gamma, Const::eps, Const::sigma, Const::T, Const::imd, Const::std_dev, Const::kb;

int main(int argc, char* argv[]){
    //command line arguments
    int n = std::stoi(argv[1]);
    double T = std::stod(argv[2]);
    double dt = std::stod(argv[3]);

    //load constants from constants.txt
    Const::load("../constants.txt");
    
    //Initialize random seed
    std::random_device rd;
    std::mt19937 gen(rd());
    
    //Dynamically allocate n particles
    Particle* particles = new Particle[n];
    for(int i = 0; i < n; i++){
        // Initialize the particles
        particles[i] = Particle(i + 1, gen);
    }
    
    size_t grid_size_x = 10;
    size_t grid_size_y = 10;
    size_t grid_size_z = 10;

   int spacing = 1.0; // Adjust the spacing as needed, works as a cell size .. we can increase this

   int box_width=10;
   int no_cells= box_width/spacing;

    // Put particles in a grid
    int particle_index = 0;
    for (size_t z = 0; z < grid_size_z; z++) {
        for (size_t y = 0; y < grid_size_y; y++) {
            for (size_t x = 0; x < grid_size_x; x++) {
                if (particle_index >= n) {
                    // All particles have been placed in the grid
                    break;
                }

                // Calculate the position of the particle in the grid
                double pos_x = x * spacing;
                double pos_y = y * spacing;
                double pos_z = z * spacing;

                // Set the position of the particle
            
                particles[particle_index].pos = Vector3d(pos_x, pos_y, pos_z);

                particle_index++;
            }
        }
    }

    extended_list* list_cells= new extended_list[no_cells * no_cells * no_cells];
    extended_list* list_particles= new extended_list[n];


    int step = 0;
    for(double t = 0; t<=T; t += dt){

        for (int i = 0; i < no_cells * no_cells * no_cells; i++) {
            list_cells[i].integer = -1;
            list_cells[i].ID = i ;
        }

        for (int i = 0; i < n; i++) {
            list_particles[i].integer = -2;
            list_particles[i].ID = i;
        }

        // std::cout<<"time:"<<t<<" ";
        // Compute force after every time step using leapfrog integration scheme
        compute_position(particles, n, dt,grid_size_x, grid_size_y, grid_size_z,list_cells,list_particles,no_cells,box_width,spacing);
        // for (size_t i = 0; i < n; i++) {
        //         std::cout << "Cell " << i << ": " << list_cells[i].integer << std::endl;
        //     }

        compute_acceleration(particles, n, dt, grid_size_x, grid_size_y, grid_size_z, spacing,list_cells, list_particles);
        // compute_acceleration(particles, n, dt,grid_size_x, grid_size_y, grid_size_z,list_cells,list_particles);
        compute_velocity(particles, n, dt);
        if(step%10 == 0)
            write_vtk(particles, n, step, grid_size_x, grid_size_y, grid_size_z);
        step++;
    }

    delete[] particles;
    delete[] list_cells;
    delete[] list_particles;
}
