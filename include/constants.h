#pragma once

#include <string>
#include <fstream>
#include <cmath>

struct Const{
    private:
        Const() {}
    public:
        static double M, R, K, gamma, eps, sigma, T, imd, kb;
        static double std_dev;
        static void load(const std::string& filename);//Load Constants from text file
};

