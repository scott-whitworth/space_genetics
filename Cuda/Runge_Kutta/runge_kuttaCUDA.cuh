#ifndef RUNGE_KUTTA_CUDA_CUH
#define RUNGE_KUTTA_CUDA_CUH

#include "../Thrust_Files/calcFourier.h"
#include "../Motion_Eqns/motion_equations.h" // Utility functions for calc_k()
#include "../Genetic_Algorithm/child.h"
#include "../Runge_Kutta/gpuMem.cuh" 

// Calls rk4SimpleCUDA to calculate the final speed and position of the rocket relative to the target
// Copies necessary values from the CPU to the GPU
// Input: calcPerS - not currently used, but is a diagnostic tool that calculates how many times the Runge Kutta algorithm ran in the kernel per second
//        generation - the children that must have the final positions calculated
//        cConstants - to access thread_block_size 
//        gpuValues - the place where all the GPU memory is allocated, stored, and freed
//                    values stored here are used in rk4SimpleCUDA calculations
//        timeInitial - the time this mission starts 
//                      this will become increasingly important when a return trip is implemented
// Output: After calling rk4SimpleCUDA, all the children in generation have their final positions and error statuses set
void callRK( double & calcPerS, Child *generation, const cudaConstants* cConstant, GPUMem & gpuValues, double timeInitial);

// the simple version of the runge_kutta algorithm, on GPU
// Finds the speed and position of an individual
// Input: children - an array of children equal to the number of threads there are that need to have their final positions calculated
//        timeInitial - the time when the spacecraft leaves its initial location (usually 0, but will not be for return trips)
//        startStepSize - the initial size of a step, which starts the same for each Child
//        absTolInput - used to calculate the error which is used to determine the next stepSize
//        n - the number of threads/ children to be processed
//        cConstant - allows you to access the thruster_type, wet_mass, and things like that
//        marsLaunchCon - the positions of Mars that have been calculated so far
// Output: The individuals in children all have their final positions and error statuses set
__global__ void rk4SimpleCUDA(Child *children, double *timeInitial, double *startStepSize, double *absTolInput, int n, const cudaConstants* cConstant, elements<double> * marsLaunchCon);


#include "runge_kuttaCUDA.cu"

#endif