#pragma once

#include <cmath>

struct Vector3d
{
    double x, y, z;

    __host__ __device__ Vector3d();

    __host__ __device__ Vector3d(double i, double j, double k) : x(i), y(j), z(k){}
 
    __host__ __device__ void set(double i, double j, double k); //set value of vector

    __host__ __device__ double norm();//Euclidean length of vector

    __host__ __device__ friend double dot(const Vector3d& v1, const Vector3d& v2);//dot product of 2 vectors

    //Operator overloading for vector calculations
    __host__ __device__ friend Vector3d operator+(const Vector3d& v1, const Vector3d& v2);

    __host__ __device__ friend Vector3d operator-(const Vector3d& v1, const Vector3d& v2);

    __host__ __device__ Vector3d operator*(const double& scalar);

    __host__ __device__ friend Vector3d operator*(const double& scalar, const Vector3d& v2);

    __host__ __device__ Vector3d operator/(const double& scalar);

    __host__ __device__ Vector3d operator>(const Vector3d& v2);

    __host__ __device__ Vector3d operator<(const Vector3d& v2);
};