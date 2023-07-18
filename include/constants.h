#pragma once

#include <string>
#include <fstream>
#include <cmath>

struct Const{
    double M, R, K, gamma, eps, sigma, 
            T, imd, kb, std_dev,gamma_t, mu;

    Const() = default;
    Const(const std::string& filename); 
};

