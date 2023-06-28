#pragma once

#include <cmath>

struct Vector3d
{
    double x, y, z;

    Vector3d();

    Vector3d(double i, double j, double k) : x(i), y(j), z(k){}
 
    void set(double i, double j, double k); //set value of vector

    double norm();//Euclidean length of vector

    friend double dot(const Vector3d& v1, const Vector3d& v2);//dot product of 2 vectors

    //Operator overloading for vector calculations
    friend Vector3d operator+(const Vector3d& v1, const Vector3d& v2);

    friend Vector3d operator-(const Vector3d& v1, const Vector3d& v2);

    Vector3d operator*(const double& scalar);

    friend Vector3d operator*(const double& scalar, const Vector3d& v2);

    Vector3d operator/(const double& scalar);
};