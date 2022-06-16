#ifndef THRUSTER_H
#define THRUSTER_H

#include "..\Config_Constants\config.h"

// file path for all data used in this file: \\fs1\phys\sankaranResearch\2019-Lauren-Mateo\Thruster

// sets starting values as given in the 2017 data excel sheet
template <class T> struct thruster {
    T P0;       // inital power in
    int type;
    T coastThreshold;

    // Constructor - takes cudaConstants to determine thruster type
    __host__ __device__ thruster<T>(const cudaConstants* gConfig);

    // //overload the stream output for elements used for writing to a file
    // template <class T> friend std::ostream & operator<<(std::ostream & Str, const thruster<T> & e); 

    // Evaluates the spacecraft's effciency for a certian iteration based off of a best fit line from file above of eff vs. power in (W).
    // Parameters:
    //         Pin: Given a certian value of power in to the spacecraft, a specific efficiency will be calculated for an iteration.
    // output: spacecraft's effciency for a certian iteration
   __host__ __device__ T calc_eff(const T & Pin);

    // Evaluates the spacecraft's fuel flow rate (mDot) for a certian iteration based off of "if statements".
    // Parameters:
    //         Pin: Given a certian value of power in to the spacecraft, a certian fuel flow rate will be set for an iteration.
    // output: spacecraft's mDot for a certian iteration
    __host__ __device__ T calc_m_Dot(const T & Pin);

    // thruster type enumeration, used in readibility for the types of thrusters rather than reading hard-coded number values
    enum THRUST_TYPE {
        NO_THRUST = 0,
        NEXT_C = 1
    };
    
    private:
        T NEXTP0 = 7330; // initial power (W)
        T NEXTm_Dot0 = 5.73E-06; // inital fuel flow rate (kg/s)

};

template <class T> __host__ __device__ T calc_accel(const T & radius, const T & z, thruster<T> & thrusterType, T & massExpelled, const T & deltaT, const bool & thrusting, const T & wetMass, const cudaConstants* cConstants);

#include "thruster.cpp"
#endif
