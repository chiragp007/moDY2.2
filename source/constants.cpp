#include "../include/constants.h"

Const::Const(const std::string& filename){
    std::fstream file;
    file.open(filename);
    std::string word;
    double var[9];
    for(int i = 1; i<=18 && file >> word; i++){
        if(i%2 == 0)
            var[i/2-1] = std::stod(word);
    }

    M = var[0];
    R = var[1];
    K = var[2];
    gamma = var[3];
    eps = var[4];
    sigma = var[5];
    T = var[6];
    imd = var[7];
    kb = var[8];
    std_dev = sqrt(kb * T/M);

    file.close();
}

