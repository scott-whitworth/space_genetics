#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H
#include "../Thrust_Files/coefficients.h"
#include "../Thrust_Files/thruster.h" // used to pass in a thruster type
#include "../Config_Constants/constants.h" // used for simulation status


// Three variations of fifth-order Runge-Kutta algorthim for system of ODEs defined in ODE45.h
// 1. rk4sys
// 2. rk4Simple
// 3. rk4Reverse

// Using adaptive time stepping  
// Parameters: 
//      timeInitial: start time (s)
//      timeFinal: end time (s)
//      times:an array that contains the time for each data point
//      y0: initial conditions (position,velocity)
//      stepSize: first time interval between data points (s)
//      y: an array which contains the soultion to the dependent variable
//      absTol: Sets the error tolerance for Runge-Kutta
//              It is very possible to set this lower than what can be represented as a double (~16 digits)
//      coeff: passes the structure containing the fourier coefficients for gamma and tau
//      accel: spacecraft's acceleration (au/s^2)
//      wetMass: mass of the spacecraft including fuel (kg)

// 1.
    // UNUSED AS OF SUMMER 2023
    // Extra parameters:
    //      Primarily used for file output
    //      lastStep: returns the index of the last element of y
    //      gamma: an array which contains all gamma values for a given run
    //      tau: an array which contains all tau values for a given run
    //      accel_output: an array which contains all accel values for a given run
    //      fuelSpent: output array to hold the aggregate fuel spent at each time step
    // Output: A dynamic array of position and velocity sets, last entry is final conditions
    //Called by trajectoryPrint()
template <class T> void rk4sys(const T & timeInitial, const T & timeFinal, T *times, const elements<T> & y0, T stepSize, elements<T> *y_new, 
                                const T & absTol, coefficients<T> coeff, T *gamma,  T *tau, int & lastStep, T *accel_output, T *fuelSpent, const T & wetMass, const cudaConstants* cConstant, elements<T> *marsIndex);

//3.
    // Comment on stepsize: Expected to be negative due to reverse integration
    // Output: writes in y the initial position of earth at the time of spacecraft launch based on an optimized trip time
    // To improve efficiency, the rk4 with single returns were split into two functions to avoid "if" statements, which are not preffered in CUDA.
template <class T> void rk4Reverse(const T & timeInitial, const T & timeFinal, const elements<T> & y0, 
                                    T stepSize, elements<T> & y_new, const T & absTol, const cudaConstants * cConstants);

// calculates k values 1 - 7 from equation and uses k values to find current and previous values of y
// error = y_new - y_prev, calculated analytically using k values
// error used in calc_scalingFactor
template <class T> __host__ __device__ void rkCalc(T & curTime, const T & timeFinal, T stepSize, elements<T> & y_new, coefficients<T> & coeff, const T & accel, 
                                                    elements<T> & error, const elements<T> & mars, T & marsCraftDist);

template <class T> void rkCalcPlanet(T & curTime, const T & timeFinal, T stepSize, elements<T> & y_new, elements<T> & error);


/**********************************************************************************************************************************/

// Calculates the scaling factor for the stepSize
// Parameters:
//      previous: The previous result of the Runge-Kutta 
//      difference: The new result minus the previous result (v-u; error from rkCalc)
//      absTol: Sets the error tolerance for Runge-Kutta
//      precThresh: passed to pmLimitCheck to set as threshold for too small
// Output: Unitless scaling coefficient which changes the time step each iteration
template <class T> __host__ __device__ T calc_scalingFactor(const elements<T> & previous , const elements<T> & difference, const T & absTol, const double precThresh);

// Error magnitude check
// Issue occurs when error is too small for a double to precicely represent
// Returns true if pmError is of reasonable accuracy
//      All components of error are within 12 orders of magnitude
// Returns false if any of the error components are more than 12 orders of magnitude from previous
// This is used to determine if the return of calc_scaling_factor should be calculated
template <class T> __host__ __device__ bool pmLimitCheck(const elements<T> & pmError, const double precThresh);

#include "runge_kutta.cpp"
#endif