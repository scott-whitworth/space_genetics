#ifndef RUNGE_KUTTA_CUDA_CUH
#define RUNGE_KUTTA_CUDA_CUH

#include "../Thrust_Files/calcFourier.h"
#include "../Motion_Eqns/motion_equations.h" // Utility functions for calc_k()
#include "../Genetic_Algorithm/child.h"
#include "../Runge_Kutta/gpuMem.cuh" 

//OUTDATED COMMENTS/CALLRK HEADER
// sets up parameters and allocates memory for and then calls rk4SimpleCUDA()
// Called by optimize() in optimization.cu
//void callRK(const int numThreads, const int blockThreads, Child *generation, double timeInitial, double stepSize, double absTol, double & calcPerS, const cudaConstants* cConstant, PlanetInfo *marsLaunchCon);

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
__global__ void rk4SimpleCUDA(Child *children, double *timeInitial, double *startStepSize, double *absTolInput, int n, const cudaConstants* cConstant, elements<double> * marsLaunchCon);


#include "runge_kuttaCUDA.cu"

#endif