#include <math.h>
#include <iostream>
#include <fstream> // for outputing to .csv file
#include <chrono>
#include <algorithm> // sort(), shuffle()
#include <random>

#include "..\Runge_Kutta\rkParameters.h"
#include "..\Config_Constants\config.h"


// Called by optimize() in optimization.cu
//void callRK(const int numThreads, const int blockThreads, Child *generation, double timeInitial, double stepSize, double absTol, double & calcPerS, const cudaConstants* cConstant, PlanetInfo *marsLaunchCon) {
void callRK(double & calcPerS, Child *generation, const cudaConstants* cConstant, GPUMem & gpuValues, double timeInitial) {

    const int numThreads = gpuValues.numThreads;
    const int blockThreads = cConstant->thread_block_size;
    double stepSize = (cConstant->triptime_min - timeInitial)/cConstant->max_numsteps;
    

    cudaEvent_t kernelStart, kernelEnd;
    cudaEventCreate(&kernelStart);
    cudaEventCreate(&kernelEnd);
    //cudaMalloc((void**) &gpuValues.devGeneration, numThreads * sizeof(Child));

    // Child *devGeneration; 
    // double *devTimeInitial;
    // double *devStepSize;
    // double *devAbsTol;
    // cudaConstants *devCConstant;
    // PlanetInfo *devMarsLaunchCon;


//     // allocate memory for the parameters passed to the device
//     cudaMalloc((void**) &devGeneration, numThreads * sizeof(Child));
//     cudaMalloc((void**) &devTimeInitial, sizeof(double));
//     cudaMalloc((void**) &devStepSize, sizeof(double)); 
//     cudaMalloc((void**) &devAbsTol, sizeof(double));
//     cudaMalloc((void**) &devCConstant, sizeof(cudaConstants));
//    // cudaMalloc((void**) &devMarsLaunchCon, sizeof());

    // copy values of parameters passed from host onto device
    cudaMemcpy(gpuValues.devGeneration, generation, numThreads * sizeof(Child), cudaMemcpyHostToDevice);
    cudaMemcpy(gpuValues.devTimeInitial, &timeInitial, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(gpuValues.devStepSize, &stepSize, sizeof(double), cudaMemcpyHostToDevice);
//     cudaMemcpy(devStepSize, &stepSize, sizeof(double), cudaMemcpyHostToDevice);
//     cudaMemcpy(devAbsTol, &absTol, sizeof(double), cudaMemcpyHostToDevice);
//     cudaMemcpy(devCConstant, cConstant, sizeof(cudaConstants), cudaMemcpyHostToDevice);
//    // cudaMemcpy(devMarsLaunchCon, marsLaunchCon, sizeof(), cudaMemcpyHostToDevice); //shallow memory copy
    
    // GPU version of rk4Simple()
    cudaEventRecord(kernelStart);
    rk4SimpleCUDA<<<(numThreads+blockThreads-1)/blockThreads,blockThreads>>>(gpuValues.devGeneration, gpuValues.devTimeInitial, gpuValues.devStepSize, gpuValues.devAbsTol, numThreads, gpuValues.devCConstant, gpuValues.devMarsLaunchCon);
    cudaEventRecord(kernelEnd);

    // copy the result of the kernel onto the host
    cudaMemcpy(generation, gpuValues.devGeneration, numThreads * sizeof(Child), cudaMemcpyDeviceToHost);


    // free memory from device
    // cudaFree(gpuValues.devGeneration);
    // cudaFree(devTimeInitial);
    // cudaFree(devStepSize);
    // cudaFree(devAbsTol);
    // cudaFree(devCConstant);
    // cudaFree(devMarsLaunchCon);

    float kernelT;
    
    cudaEventSynchronize(kernelEnd);

    cudaEventElapsedTime(&kernelT, kernelStart, kernelEnd);
    
    calcPerS = numThreads / (kernelT / 1000.0); // how many times the Runge Kutta algorithm ran in the kernel per second
}

// seperate conditions are passed for each thread, but timeInitial, stepSize, and absTol are the same for every thread
__global__ void rk4SimpleCUDA(Child *children, double *timeInitial, double *startStepSize, double *absTolInput, int n, const cudaConstants* cConstant, elements<double>* marsLaunchCon) {

     // pre-calculate a table of Mars' position within possible mission time range
        // defined as global variable
        // accessed on the CPU when individuals are initialized
    

    int threadId = threadIdx.x + blockIdx.x * blockDim.x;
    if (threadId < n) {
        rkParameters<double> threadRKParameters = children[threadId].startParams; // get the parameters for this thread

        elements<double> curPos = threadRKParameters.y0; // start with the initial conditions of the spacecraft
        //position of mars at the current time
        elements<double> mars;

        // storing copies of the input values
        double stepSize = *startStepSize;
        double absTol = *absTolInput;
        double curTime = *timeInitial;
        double startTime = *timeInitial;
        double curAccel = 0;

        thruster<double> thrust(cConstant);

        double massFuelSpent = 0; // mass of total fuel expended (kg) starts at 0

        bool coast; // to hold the result from calc_coast()

        elements<double> error; // holds output of previous value from rkCalc

        while (curTime < threadRKParameters.tripTime) {

            // Check the thruster type before performing calculations
            if (cConstant->thruster_type == thruster<double>::NO_THRUST) {
                coast = curAccel = 0;
            }
            else {
                coast = calc_coast(threadRKParameters.coeff, curTime, threadRKParameters.tripTime, thrust);
                curAccel = calc_accel(curPos.r, curPos.z, thrust, massFuelSpent, stepSize, coast, static_cast<double>(cConstant->wet_mass), cConstant);
            }

            //Needs to be triptime - curtime to get the correct index for mars
            //when curtime = triptime, this will give us the final position of mars at impact
            //this is because getConditionDev takes in seconds before the spacecraft reaches the tarfet
            mars = getConditionDev(threadRKParameters.tripTime - curTime, cConstant, marsLaunchCon);

            //calculate the distance between mars and the spacecraft (|R|^2)
            double marsCraftDist = sqrt(pow(mars.r, 2) + pow(curPos.r, 2) + pow(curPos.z - mars.z, 2) - (2*curPos.r*mars.r*cos(mars.theta-curPos.theta)));
            
            // calculate k values and get new value of y
            rkCalc(curTime, threadRKParameters.tripTime, stepSize, curPos, threadRKParameters.coeff, curAccel, error, mars, marsCraftDist); 

            curTime += stepSize; // update the current time in the simulation
            
            stepSize *= calc_scalingFactor(curPos-error,error,absTol, cConstant->doublePrecThresh); // Alter the step size for the next iteration

            // The step size cannot exceed the total time divided by 2 and cannot be smaller than the total time divided by 1000
            if (stepSize > (threadRKParameters.tripTime - startTime) / cConstant->min_numsteps) {
                stepSize = (threadRKParameters.tripTime - startTime) / cConstant->min_numsteps;
            }
            else if (stepSize < (threadRKParameters.tripTime - startTime) / cConstant->max_numsteps){
                stepSize = (threadRKParameters.tripTime - startTime) / cConstant->max_numsteps;
            }

            //if the spacecraft is within Mars' sphere of influence, make the stepSize as small as possible
            //  this will increase the resolution of the spacecraft's position in this time frame
            if (marsCraftDist < MSOI*cConstant->MSOI_error){
                stepSize = (threadRKParameters.tripTime - startTime) / (cConstant->MSOI_steps*cConstant->max_numsteps);
            }
            //count the steps taken for this threads calculations
            children[threadId].stepCount++;

            if ( (curTime + stepSize) > threadRKParameters.tripTime) {
                stepSize = (threadRKParameters.tripTime - curTime); // shorten the last step to end exactly at time final
            }

            // if the spacecraft is within 0.5 au of the sun, the radial position of the spacecraft artificially increases to 1000, to force that path to not be used in the optimization.
            if ( sqrt(pow(curPos.r,2) + pow(curPos.z,2)) < cConstant->sun_r_min) { //maybe issue is with using pow? I doubt it, but we could always try curPos.r*curPos.r + curPos.z*curPos.z < sun_r_min*sun_r_min?
                //This is a bad result, needs to be set to be removed
                //Setting the child's status to be a sun error
                children[threadId].errorStatus = SUN_ERROR;//Are all the children's errorStatus set to SUN_ERROR?

                return;
            }
            
            //if(sqrt(pow(mars.r, 2) + pow(curPos.r, 2) + pow(curPos.z - mars.z, 2) - (2*curPos.r*mars.r*cos(mars.theta-curPos.theta))) < cConstant->gravAssistDist){
            if(marsCraftDist < cConstant->gravAssistDist){
                //This is a bad result, needs to be set to be removed
                //Setting the child's status to be a sun error
                children[threadId].errorStatus = MARS_ERROR;//Are all the children's errorStatus set to SUN_ERROR?

                return;
            }
            
        }
        //Give the child its final calculated position
        children[threadId].finalPos = curPos;

        //if it is not a SUN_ERROR then it is valid
        children[threadId].errorStatus = VALID;

        return;
    }
    return;
}