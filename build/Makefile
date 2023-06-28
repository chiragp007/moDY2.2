FLAGS = -g -O3 -Wall -pedantic -std=c++17 

md: main.o particle.o vec.o constants.o 
	g++ $(FLAGS) -o md main.o particle.o vec.o constants.o 

main.o: ../source/main.cpp ../include/vec.h ../include/constants.h ../include/particle.h ../include/vtkwrite.h
	g++ $(FLAGS) -c ../source/main.cpp

particle.o: ../source/particle.cpp ../include/particle.h ../include/vec.h ../include/constants.h
	g++ $(FLAGS) -c ../source/particle.cpp

vec.o: ../source/vec.cpp ../include/vec.h
	g++ $(FLAGS) -c ../source/vec.cpp

constants.o: ../source/constants.cpp ../include/constants.h
	g++ $(FLAGS) -c ../source/constants.cpp

clean:
	rm -rf *.o md ../output/*.vtk


