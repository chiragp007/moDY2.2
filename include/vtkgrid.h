#pragma once

#include<iostream>
#include<fstream>

using namespace std;

void create_wireframe(
    const double xmin, const double xmax,
    const double ymin, const double ymax,
    const double zmin, const double zmax
){
    fstream f;
    f.open("../output/wireframe.vtk", ios::out);
    if(f.is_open()){

        f << "# vtk DataFile Version 2.0\n" << "Wireframe Boundry\n"
          << "ASCII\n" << "DATASET POLYDATA\n" << "POINTS 8 float\n";

        f << xmin << " " << ymin << " " << zmin << "\n";
        f << xmax << " " << ymin << " " << zmin << "\n";
        f << xmax << " " << ymax << " " << zmin << "\n";
        f << xmin << " " << ymax << " " << zmin << "\n";
        f << xmin << " " << ymin << " " << zmax << "\n";
        f << xmax << " " << ymin << " " << zmax << "\n";
        f << xmax << " " << ymax << " " << zmax << "\n";
        f << xmin << " " << ymax << " " << zmax << "\n";

        f << "LINES " << 12 << " " << 36 << endl;

        f << 2 << " " << 0 << " " << 1 << endl;
        f << 2 << " " << 1 << " " << 2 << endl;
        f << 2 << " " << 2 << " " << 3 << endl;
        f << 2 << " " << 3 << " " << 0 << endl;
        f << 2 << " " << 4 << " " << 5 << endl;
        f << 2 << " " << 5 << " " << 6 << endl;
        f << 2 << " " << 6 << " " << 7 << endl;
        f << 2 << " " << 7 << " " << 4 << endl;
        f << 2 << " " << 0 << " " << 4 << endl;
        f << 2 << " " << 1 << " " << 5 << endl;
        f << 2 << " " << 2 << " " << 6 << endl;
        f << 2 << " " << 3 << " " << 7 << endl;

        f.close();

        cout << "\033[1;32m\n\tThe file wireframe.vtk has been created.\n\033[0m\n";

    }
    else{
        cout << "\033[1;31m\n\tThe file wireframe.vtk cannot be created.\n\033[0m\n";
        exit(-1);
    }
}