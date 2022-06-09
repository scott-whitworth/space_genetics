#ifndef RUNGE_KUTTA_CUDA_CUH
#define RUNGE_KUTTA_CUDA_CUH

#include "../Thrust_FIles/calcFourier.h"
#include "../Motion_Eqns/motion_equations.h" // Utility functions for calc_k()
#include "../Genetic_Algorithm/child.h"

// sets up parameters and allocates memory for and then calls rk4SimpleCUDA()
// Called by optimize() in optimization.cu
void callRK(const int numThreads, const int blockThreads, Child *generation, double timeInitial, double stepSize, double absTol, double & calcPerS, const cudaConstants* cConstant);

// the simple version of the runge_kutta algorithm, on GPU
__global__ void rk4SimpleCUDA(Child *children, double *timeInitial, double *startStepSize, double *absTolInput, int n, const cudaConstants* cConstant);


#include "runge_kuttaCUDA.cu"
#endif