#ifndef CALCFOURIER_H
#define CALCFOURIER_H

// General fourier series function which can be called for specific parameters (gamma,tau,and coast)
// Parameters:
//         series[]: input array, the coefficients of the series
//         series_size: size of the input array
//         curTime: current time (s) for calculated gamma
//         timeFinal: end time (s), used to normalize t
// output: the value of the fourier series evaluated at a curTime
template <class T> __host__ __device__ T calc_Series(const T series[], const int series_size, const T & curTime, const T & timeFinal);

// Calculates gamma (in-plane angle) at a specific time using Fourier series coefficients
// Parameters:
//         coeff: coefficients structure, specifically the gamma components
//                coeff.gamma is an array of optimized <T> values 
//                coeff.gammaSize is the size of this array 
//         curTime: current time (s) for calculated gamma
//         timeFinal: end time (s), used to normalize t
// output: in-plane angle derived from normalized time and gamma Fourier series
template <class T> __host__ __device__ T calc_gamma(const coefficients<T> & coeff,const T & curTime, const T & timeFinal);


// Calculates tau (out-of-plane angle) at a specific time using Fourier series coefficients
// Parameters:
//         coeff: coefficients structure, specifically the tau components
//                coeff.tau is an array of optimized <T> values 
//                coeff.tauSize is the size of this array
//         curTime: current time (s) for calculated gamma
//         timeFinal: end time (s), used to normalize t
// output: in-plane angle derived from normalized time and tau Fourier series
template <class T> __host__ __device__ T calc_tau(const coefficients<T> & coeff,const T & curTime, const T & timeFinal);

// Evaluates whether the spacecraft is accelerating or coasting for a specific iteration
// Parameters:
//         coeff: coefficients structure, specifically the tau components
//                coeff.coast is an array of optimized <T> values 
//                coeff.coastSize is the size of this array
//                coeff.coastThreshold is the value at which above, the spacecraft accelerates and below, the spacecraft coasts
//         curTime: current time (s) for calculated gamma
//         timeFinal: end time (s), used to normalize t
// output: 
//         above the threshold: returns thrusting = 1
//         below the threshold: returns thrusting = 0
template <class T> __host__ __device__ bool calc_coast(coefficients<T> & coeff, const T & curTime, const T & timeFinal, thruster<T> & thrust);

#include "calcFourier.cpp"
#endif