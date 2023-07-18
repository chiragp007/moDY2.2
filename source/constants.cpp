#include <sstream> // Add this include for stringstream operations

#include "../include/constants.h"

Const::Const(const std::string& filename) {
    std::fstream file;
    file.open(filename);

    // Read all constants using string stream
    std::string line;
    int i = 0;
    while (std::getline(file, line) && i < 11) {
        std::stringstream ss(line);
        std::string key;
        ss >> key;
        if (key == "M")
            ss >> M;
        else if (key == "R")
            ss >> R;
        else if (key == "K")
            ss >> K;
        else if (key == "gamma")
            ss >> gamma;
        else if (key == "eps")
            ss >> eps;
        else if (key == "sigma")
            ss >> sigma;
        else if (key == "T")
            ss >> T;
        else if (key == "imd")
            ss >> imd;
        else if (key == "kb")
            ss >> kb;
        else if (key == "gamma_r")
            ss >> gamma_t;
        else if (key == "mu")
            ss >> mu;

        i++;
    }
    std_dev = sqrt(kb * T / M);

    file.close();
}
