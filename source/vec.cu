#include "../include/vec.cuh"

__host__ __device__ Vector3d::Vector3d() {
    x = 0; y = 0; z = 0;
}

__host__ __device__ void Vector3d::set(double i, double j, double k){
    x = i; y = j; z = k;
}

__host__ __device__ Vector3d operator+(const Vector3d& v1, const Vector3d& v2){
    Vector3d v;
    v.set(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
    return v;
}

__host__ __device__ Vector3d operator-(const Vector3d& v1, const Vector3d& v2){
    Vector3d v;
    v.set(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
    return v;
}

__host__ __device__ Vector3d cross(const Vector3d& v1, const Vector3d& v2){
    return Vector3d(v1.y * v2.z + v1.z * v2.y, -(v1.x * v2.z + v1.z * v2.x), v1.x * v2.y + v1.y * v2.x);
}

__host__ __device__ Vector3d Vector3d::operator*(const double& scalar){
    Vector3d v;
    v.set(scalar * this->x, scalar * this->y, scalar * this->z);
    return v;
}

__host__ __device__ Vector3d operator*(const double& scalar, const Vector3d& v2){
    Vector3d v;
    v.set(scalar * v2.x, scalar * v2.y, scalar * v2.z);
    return v;
}

__host__ __device__ Vector3d Vector3d::operator/(const double& scalar){
    Vector3d v;
    v.set(this->x / scalar, this->y / scalar, this->z / scalar);
    return v;
}

__host__ __device__ double Vector3d::norm(){
    return sqrt(pow(this->x, 2.0) + pow(this->y, 2.0) + pow(this->z, 2.0));
}

__host__ __device__ double dot(const Vector3d& v1, const Vector3d& v2){
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

__host__ __device__ Vector3d Vector3d::operator>(const Vector3d& v2){
    return Vector3d(this->x > v2.x, this->y > v2.y, this->z > v2.z);
}

__host__ __device__ Vector3d Vector3d::operator<(const Vector3d& v2){
    return Vector3d(this->x < v2.x, this->y < v2.y, this->z < v2.z);
}