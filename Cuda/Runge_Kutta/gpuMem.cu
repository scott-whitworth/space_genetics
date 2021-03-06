#include "gpuMem.cuh"

 void GPUMem::initialize(const cudaConstants* cConstants, const int & marsConSize, const elements<double>* marsLaunchCon){
    //all these initializations taken from (old) line 31 of genetic_algorithm.cpp in callRK
    //stepSize = (cConstants->orbitalPeriod / cConstants->max_numsteps);
    absTol = cConstants->rk_tol;
    numThreads = cConstants->num_individuals;
    
    // allocate memory for the parameters passed to the device
    cudaMalloc((void**) &devGeneration, numThreads * sizeof(Child)); 
    cudaMalloc((void**) &devTimeInitial, sizeof(double));
    cudaMalloc((void**) &devStepSize, sizeof(double)); 
    cudaMalloc((void**) &devAbsTol, sizeof(double));
    cudaMalloc((void**) &devCConstant, sizeof(cudaConstants));
    cudaMalloc((void**) &devMarsLaunchCon, marsConSize * sizeof(double));

    // copy values of parameters passed from host onto device
    //cudaMemcpy(devStepSize, &stepSize, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devAbsTol, &absTol, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devCConstant, cConstants, sizeof(cudaConstants), cudaMemcpyHostToDevice);
    cudaMemcpy(devMarsLaunchCon, marsLaunchCon, marsConSize * sizeof(double), cudaMemcpyHostToDevice);
}

void GPUMem::free(){
    cudaFree(devGeneration);
    cudaFree(devTimeInitial);
    cudaFree(devStepSize);
    cudaFree(devAbsTol);
    cudaFree(devCConstant);
    cudaFree(devMarsLaunchCon);
}