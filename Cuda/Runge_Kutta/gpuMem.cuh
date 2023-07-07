#ifndef GPU_MEM_CUH
#define GPU_MEM_CUH

#include "../Genetic_Algorithm/child.h"
#include "..\Runge_Kutta\rkParameters.h"
#include "..\Config_Constants\config.h"

//This struct has all the device parameters need for the callRK process
struct GPUMem
{
    // double stepSize; 
    double absTol;
    int numThreads;

    Child *devGeneration; 
    // double *devTimeInitial;
    // double *devStepSize;
    double *devAbsTol;
    cudaConstants *devCConstant;
    elements<double> *devMarsLaunchCon;

    double *devTime_steps;
    elements<double> *devY_steps;
    double *devGamma_steps;
    double *devTau_steps;
    double *devAccel_steps;
    double *devFuel_steps;

    //initializes the device parameters and allocates the space on the GPU
    //input: cConstants - needed to allocate space
    //       marsConSize - size of (marsCon * elements<double>)
    //       marsLaunchCon - iformation regarding mars' position to be set on GPU
    //output: space allocated on GPU and information stored in that space for callRK
    void initialize(const cudaConstants* cConstants, const int & marsConSize, const elements<double>* marsLaunchCon);

    //deallocates memory stored to GPU
    void free();

};


#include "gpuMem.cu"

#endif