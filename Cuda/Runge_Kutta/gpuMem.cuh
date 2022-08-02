#ifndef GPU_MEM_CUH
#define GPU_MEM_CUH

#include "../Genetic_Algorithm/child.h"
#include "..\Runge_Kutta\rkParameters.h"
#include "..\Config_Constants\config.h"

//This struct has all the device parameters need for the callRK process
struct GPUMem
{
    double stepSize; 
    double absTol;
    int numThreads;

    Child *devGeneration; 
    double *devTimeInitial;
    double *devStepSize;
    double *devAbsTol;
    cudaConstants *devCConstant;
    elements<double> *devMarsLaunchCon;
    
    void initialize(const cudaConstants* cConstants, const int & marsConSize, const elements<double>* marsLaunchCon);

    void free();

};


#include "gpuMem.cu"

#endif