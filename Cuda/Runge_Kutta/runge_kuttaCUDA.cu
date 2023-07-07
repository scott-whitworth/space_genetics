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

    //Counters for the number of children who are finished being simulated and how many simulations the code has ran for this generation
    int completedNum;
    int simulationNum = 0;

    const int numThreads = gpuValues.numThreads;
    const int blockThreads = cConstant->thread_block_size;
    double stepSize = (cConstant->triptime_min - timeInitial)/cConstant->max_numsteps;
    
    //Do while loop will keep simulating while there are still children to simulate
    do {

        //Reset the count of children who have finished their simulation
        completedNum = 0;

        cudaEvent_t kernelStart, kernelEnd;
        cudaEventCreate(&kernelStart);
        cudaEventCreate(&kernelEnd);

        // copy values of parameters passed from host onto device
        cudaMemcpy(gpuValues.devGeneration, generation, numThreads * sizeof(Child), cudaMemcpyHostToDevice);
        // cudaMemcpy(gpuValues.devTimeInitial, &timeInitial, sizeof(double), cudaMemcpyHostToDevice);
        // cudaMemcpy(gpuValues.devStepSize, &stepSize, sizeof(double), cudaMemcpyHostToDevice);

        // GPU version of rk4Simple()
        cudaEventRecord(kernelStart);
        rk4CUDASim<<<(numThreads+blockThreads-1)/blockThreads,blockThreads>>>(gpuValues.devGeneration, gpuValues.devAbsTol, numThreads, gpuValues.devCConstant, gpuValues.devMarsLaunchCon, gpuValues.devTime_steps, gpuValues.devY_steps, gpuValues.devGamma_steps, gpuValues.devTau_steps, gpuValues.devAccel_steps, gpuValues.devFuel_steps);
        cudaEventRecord(kernelEnd);

        // copy the result of the kernel onto the host
        cudaMemcpy(generation, gpuValues.devGeneration, numThreads * sizeof(Child), cudaMemcpyDeviceToHost);

        //Output how many simulation cycles the code has done
        //std::cout << "\nCompleted simulation cycle " << simulationNum << " for this generation.\n";
        std::cout << "_";

        simulationNum++; //Add one to simulation num for the next cycle

        float kernelT;
        
        cudaEventSynchronize(kernelEnd);

        cudaEventElapsedTime(&kernelT, kernelStart, kernelEnd);
        
        calcPerS = numThreads / (kernelT / 1000.0); // how many times the Runge Kutta algorithm ran in the kernel per second

        //Go through simulated individuals to see how many are finished
        for (int i = 0; i < cConstant->num_individuals; i++)
        {
            if(generation[i].simStatus == COMPLETED_SIM) {
                completedNum++;
            }
        }
    } while (completedNum < cConstant->num_individuals);
}

// seperate conditions are passed for each thread, but timeInitial, stepSize, and absTol are the same for every thread
__global__ void rk4CUDASim(Child *children, double *absTolInput, int n, const cudaConstants* cConstant, elements<double> * marsLaunchCon, double *time_steps, elements<double> *y_steps, double *gamma_steps, double *tau_steps, double *accel_steps, double *fuel_steps) {
   
   int threadId = threadIdx.x + blockIdx.x * blockDim.x;
    if (threadId < n) {
        // children[threadId].simNum++;
        //Check if this child has been simulated already
        if (children[threadId].simStatus != COMPLETED_SIM) {
            //If not, run it's simulation

            //Check to see if the child is about to be simulated too many times
            if ((children[threadId].simNum + 1) > cConstant->maxSimNum) {
                //If not, run the next simulation and increase the sim num to reflect this
                children[threadId].simNum++;

                //Assign an error to the child because it has been running for too long
                children[threadId].errorStatus = SIMNUM_ERROR;

                //Mark it as having completed it's run
                children[threadId].simStatus = COMPLETED_SIM;

                //Quit the simulation
                return;
            }
            
            rkParameters<double> threadRKParameters = children[threadId].startParams; // get the parameters for this thread

            // storing copies of the input values
            double stepSize;
            double absTol = *absTolInput;
            double startTime;
            double endTime;
            double curAccel = 0;

            //Stores the total mass of the fuel spent during the simulation
            double massFuelSpent;

            //Initial time and position of the child
            double curTime;
            elements<double> curPos;

            //Stores the angular momentum if the individual is entering an SOI
            //  Used to check if the assist was bad
            double soiEntryh;

            //Set the initial curTime and curPos depending on if the child has been ran
            if (children[threadId].simStatus == INITIAL_SIM) {
                //Child has not been simulated, set the initial curTime to the start time of the simulation
                //Set at 0 initially as no time has passed within the simulation yet
                curTime = 0;
                //Set the start time to the total trip time
                startTime = curTime;
                // start with the initial conditions of the spacecraft
                curPos = threadRKParameters.y0; 

                massFuelSpent = 0; // mass of total fuel expended (kg) starts at 0
            }
            else {
                //Child has been partially simulated (means it has entered a SOI), set initial curTime to the child's simStartTime variable
                //That will be set to the simStartTime variable for child
                curTime = children[threadId].simStartTime;
                //Set the start time to this simulation's start time
                startTime = children[threadId].simStartTime;
                //Get the child's simStartPos, which will have the elements of the child at the last step of the last simulation
                curPos = children[threadId].simStartPos;

                //Get the mass of how much fuel the child has spent on previous simulations
                massFuelSpent = children[threadId].fuelSpent;
            }

            //Check to see if this simulation occurs within a sphere of influence
            //This determines the endTime
            if(children[threadId].simStatus == INSIDE_SOI) {
                //If this is a mission inside a SOI, set the final to be a gravAssistTime difference from the start time
                endTime = startTime + cConstant->gravAssistTime;
                //endTime = startTime + ((threadRKParameters.tripTime - startTime) * cConstant->gravAssistTimeFrac);

                //Calculate the entry angular momentum
                soiEntryh = threadRKParameters.y0.r * threadRKParameters.y0.vtheta;

                //Check to make sure that endTime is not further than tripTime
                if(endTime > threadRKParameters.tripTime) {
                    endTime = threadRKParameters.tripTime;
                }
            }
            //In all other scenerios, the end time is the triptime
            else {
                endTime = threadRKParameters.tripTime;
            }

            thruster<double> thrust(cConstant);

            bool coast; //=1 means thrusting, from calc_coast()

            elements<double> error; // holds output of previous value from rkCalc

            stepSize = (endTime - startTime) / cConstant->max_numsteps;

            while (curTime < endTime) {
                time_steps[children[threadId].stepCount] = curTime;

                // Check the thruster type before performing calculations
                if (cConstant->thruster_type == thruster<double>::NO_THRUST) {
                    coast = curAccel = 0;
                }
                else {
                    coast = calc_coast(threadRKParameters.coeff, curTime, threadRKParameters.tripTime, thrust);

                    curAccel = calc_accel(curPos.r, curPos.z, thrust, massFuelSpent, stepSize, coast, static_cast<double>(cConstant->wet_mass), cConstant);

                    gamma_steps[children[threadId].stepCount] = calc_gamma(children[threadId].startParams.coeff, curTime, threadRKParameters.tripTime);

                    tau_steps[children[threadId].stepCount] = calc_tau(children[threadId].startParams.coeff, curTime, threadRKParameters.tripTime);

                    children[threadId].fuelSpent = massFuelSpent;
                }

                accel_steps[children[threadId].stepCount] = curAccel;
                fuel_steps[children[threadId].stepCount] = massFuelSpent;

                //Needs to be triptime - curtime to get the correct index for mars
                //when curtime = triptime, this will give us the final position of mars at impact
                //this is because getConditionDev takes in seconds before the spacecraft reaches the target
                elements<double> mars = getConditionDev(threadRKParameters.tripTime - curTime, cConstant, marsLaunchCon);

                //calculate the distance between mars and the spacecraft (|R|^2)
                double marsCraftDist = sqrt(pow(mars.r, 2) + pow(curPos.r, 2) + pow(curPos.z - mars.z, 2) - (2*curPos.r*mars.r*cos(mars.theta-curPos.theta)));

                //See if the child is closest it has been to Mars so far this run
                //This is only updated if Mars is in betweem the craft amd the target
                if (marsCraftDist < children[threadId].minMarsDist) {
                    children[threadId].minMarsDist = marsCraftDist;
                }
                
                // calculate k values and get new value of y
                rkCalc(curTime, threadRKParameters.tripTime, stepSize, curPos, threadRKParameters.coeff, curAccel, error, mars, marsCraftDist); 

                curTime += stepSize; // update the current time in the simulation
                
                //stepSize *= calc_scalingFactor(curPos-error,error,absTol, cConstant->doublePrecThresh); // Alter the step size for the next iteration

                //// The step size cannot exceed the total time divided by 2 and cannot be smaller than the total time divided by 1000
                //if (stepSize > (endTime - startTime) / cConstant->min_numsteps) {
                //    stepSize = (endTime - startTime) / cConstant->min_numsteps;
                //}
                //else if (stepSize < (endTime - startTime) / cConstant->max_numsteps){
                //    stepSize = (endTime - startTime) / cConstant->max_numsteps;
                //}

                if ( (curTime + stepSize) > endTime) {
                    stepSize = (endTime - curTime); // shorten the last step to end exactly at time final
                }

                // if the spacecraft is within 0.5 au of the sun, the radial position of the spacecraft artificially increases to 1000, to force that path to not be used in the optimization.
                if ( sqrt(pow(curPos.r,2) + pow(curPos.z,2)) < cConstant->sun_r_min) { //maybe issue is with using pow? I doubt it, but we could always try curPos.r*curPos.r + curPos.z*curPos.z < sun_r_min*sun_r_min?
                    //This is a bad result, needs to be set to be removed
                    //Setting the child's status to be a sun error
                    children[threadId].errorStatus = SUN_ERROR;//Are all the children's errorStatus set to SUN_ERROR?

                    //Set the child's sim status to complete so it ins't ran any further
                    children[threadId].simStatus = COMPLETED_SIM; 

                    return;
                }
                
                //Check to see if the child is too close to Mars
                if(marsCraftDist < cConstant->gravAssistDist){
                    //This is a bad result, needs to be set to be removed
                    //Setting the child's status to be a mars error
                    children[threadId].errorStatus = MARS_ERROR;
                    
                    //Set the child's sim status to complete so it ins't ran any further
                    children[threadId].simStatus = COMPLETED_SIM;

                    return;
                }

                //Check to see if the curTime is less than triptime
                if (curTime < endTime) {
                    //If so, check to see if the child triggers the conditions for a new run

                    //See if child has entered MSOI
                    if (marsCraftDist < MSOI*cConstant->MSOI_scale && children[threadId].simStatus != INSIDE_SOI) {
                        //Set the final time and position for this portion of the child's simulation
                        children[threadId].simStartTime = curTime;
                        children[threadId].simStartPos = curPos;

                        //Set the child status to inside SOI
                        children[threadId].simStatus = INSIDE_SOI;

                        return;
                    }

                    //Check if child has exited MSOI after being inside it
                    if (marsCraftDist > MSOI*cConstant->MSOI_scale && children[threadId].simStatus == INSIDE_SOI) {
                        //Set the final time and position for this portion of the child's simulation
                        children[threadId].simStartTime = curTime;
                        children[threadId].simStartPos = curPos;

                        //Set the child status to outide SOI
                        children[threadId].simStatus = OUTSIDE_SOI;

                        children[threadId].orbithChange = (curPos.r * curPos.vtheta) - soiEntryh;

                        //Check to make sure the orbithChange isn't less than 0 
                        if (children[threadId].orbithChange < 1e-14) {
                            children[threadId].orbithChange = 5e-15;
                        }

                        //Check to see if this was a bad assist
                        // if (curPos.r * curPos.vtheta < soiEntryh) {
                        //     //If it is a bad assist, set the error
                        //     children[threadId].errorStatus = BAD_ASSIST;
                        //     children[threadId].simStatus = COMPLETED_SIM;
                        // }

                        return;
                    }

                    y_steps[children[threadId].stepCount] = curPos;
                
                    //count the steps taken for this threads calculations
                    children[threadId].stepCount++;

                }  //end if (curTime < endTime)
            } // end while (curTime < endTime)

            //Check to see if this simulation has completed its total runtime
            if (endTime >= threadRKParameters.tripTime) {
                //If here, the individual's simulation has completely ended, regardless of its current simStatus
                
                //Give the child its final calculated position
                children[threadId].finalPos = curPos;

                //if it is not a SUN_ERROR then it is valid
                children[threadId].errorStatus = VALID;

                //The child has finished it's run, set the simStatus to completed
                children[threadId].simStatus = COMPLETED_SIM;

                //children[threadId].fuelSpent = massFuelSpent;
            }
            //else: endTime != tripTime
            else {
                //This means the simulation was for an individual in a SOI, but it never escaped the SOI in the estimated time
                //The individual will need be run through a SOI simulation again 
                //The sim status would already be INSIDE_SOI, so only the position and time need to be recorded
                children[threadId].simStartPos = curPos;
                children[threadId].simStartTime = curTime;  
            }

            return;
        }
    }
}