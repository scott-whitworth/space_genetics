#ifndef COEFFICIENTS_H
#define COEFFICIENTS_H

#include "../Config_Constants/config.h"
#include <iostream>

//Utility to calculate the coefficients vector
//Coefficients struct holds gamma and tau values
template <class T> struct coefficients {
    template <class T> void initCoefficient(double x[], coefficients<T> & coeff, const cudaConstants* cConstants);
    
    // fourth order fourier series
    // setting the size of gamma array (in-plane coefficients angle)
    const static int gammaSize = GAMMA_ARRAY_SIZE; 
    T gamma[gammaSize]; 

    // first order fourier series
    // setting the size of tau array (out-of-plane coefficients angle)
    const static int tauSize = TAU_ARRAY_SIZE; 
    T tau[tauSize]; 

    // second order fourier series
    // setting the size of coast array
    const static int coastSize = COAST_ARRAY_SIZE;
    T coast[coastSize];

    // threshold: after sine squared of the fourier series is evaluated, if above the threshold, acceleration occurs. When below, coasting occurs.
    // Is assigned a value from cudaConstants in initCoefficient()
    T coastThreshold;

    //overload the stream output for elements used for writing to a file
    template <class U> friend std::ostream & operator<<(std::ostream & Str, const coefficients<T> & e);
};

#include "coefficients.cpp"
#endif