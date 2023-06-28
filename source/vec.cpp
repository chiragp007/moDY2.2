#include "../include/vec.h"

Vector3d::Vector3d() {
    x = 0; y = 0; z = 0;
}

void Vector3d::set(double i, double j, double k){
    x = i; y = j; z = k;
}

Vector3d operator+(const Vector3d& v1, const Vector3d& v2){
    Vector3d v;
    v.set(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
    return v;
}

Vector3d operator-(const Vector3d& v1, const Vector3d& v2){
    Vector3d v;
    v.set(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
    return v;
}

Vector3d Vector3d::operator*(const double& scalar){
    Vector3d v;
    v.set(scalar * this->x, scalar * this->y, scalar * this->z);
    return v;
}

Vector3d operator*(const double& scalar, const Vector3d& v2){
    Vector3d v;
    v.set(scalar * v2.x, scalar * v2.y, scalar * v2.z);
    return v;
}

Vector3d Vector3d::operator/(const double& scalar){
    Vector3d v;
    v.set(this->x / scalar, this->y / scalar, this->z / scalar);
    return v;
}

double Vector3d::norm(){
    return sqrt(pow(this->x, 2.0) + pow(this->y, 2.0) + pow(this->z, 2.0));
}

double dot(const Vector3d& v1, const Vector3d& v2){
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}