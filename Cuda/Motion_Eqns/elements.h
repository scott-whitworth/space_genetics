#ifndef ELEMENTS_H
#define ELEMENTS_H

#include <iostream> // used for overloading << operator

//elements struct holds k values / dependent variable values in rk4sys
template <class T> struct elements {
    // All in relation to the plane of the sun in cylindrical coordinates

    // Positions
    T r;     //radius (in plane) - AU
    T theta; //angular position (in plane) - Radians
    T z;     //axial position (out-of-plane) - AU

    // Velocities
    T vr;     // Radial velocity (in plane) - AU/s
    T vtheta; // tangential velocity (in plane) - AU/s
    T vz;     // axial velocity (out-of-plane) - AU/s

    // Constructor which sets all the components
    __host__ __device__ elements<T>(T r0, T theta0, T z0, T vr0, T vtheta0, T vz0);

    // Default constructor which takes no values and sets everything to 0
    __host__ __device__ elements<T>();

    // Overload operators to do math on all the elements in the struct seperately
    // Treating each element as a matrix operation

    // Operator overloads, used in runge kutta for the math between elements
    __host__ __device__ elements<T> operator+(const elements<T>& e) const;
    __host__ __device__ elements<T> operator-(const elements<T>& e) const;
    // __host__ __device__ elements<T> operator*(const elements<T>& e) const;
    // __host__ __device__ elements<T> operator/(const elements<T>& e) const;

    // Operator overloads, for scalars
    __host__ __device__ elements<T> operator*(const T& i) const;
    __host__ __device__ elements<T> operator/(const T& i) const;


    // Comparison method to compare this elements to another, returning true (equivalent) within a threshold
    // input: other - another elements object to compare to
    //        comp_Thresh - a threshold for comparison, so when elements are considered equivalent it is a not necessarrily 100% (when not equal to 0)
    // output: returns true if for all properties (pos & vel) the difference between them is less than or equal to comp_Thresh
    bool compare(const elements<T> & other, T comp_Thresh);

    // Overload the stream output for elements used for writing to a file
    template <class U> friend std::ostream & operator<<(std::ostream & Str, const elements<T> & e); 
};

#include "elements.cpp"

#endif