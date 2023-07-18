CC_FLAGS = -O3 -pedantic -std=c++17 
NVCC_FLAGS = -dc -lcudadevrt

CUDA_ROOT_DIR=/usr/local/cuda
CUDA_LIB_DIR= -L$(CUDA_ROOT_DIR)/lib64
CUDA_INC_DIR= -I$(CUDA_ROOT_DIR)/include
CUDA_LINK_LIBS = -lcudart

md: main.o particle.o vec.o constants.o
	nvcc main.o particle.o vec.o constants.o -o md $(CUDA_LINK_LIBS) $(CUDA_LIB_DIR) $(CUDA_LIB_DIR)

main.o: ../source/main.cu ../include/vec.cuh ../include/constants.h ../include/particle.cuh ../include/vtkwrite.h ../include/vtkgrid.h
	nvcc $(NVCC_FLAGS) -c ../source/main.cu

particle.o: ../source/particle.cu ../include/particle.cuh ../include/vec.cuh ../include/constants.h
	nvcc $(NVCC_FLAGS) -c ../source/particle.cu

vec.o: ../source/vec.cu ../include/vec.cuh
	nvcc $(NVCC_FLAGS) -c ../source/vec.cu

constants.o: ../source/constants.cpp ../include/constants.h
	g++ $(CC_FLAGS) -c ../source/constants.cpp

clean:
	rm -rf *.o *.obj *.exe  md ../output/*.vtk


