#include <string>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <direct.h>
//#include <algorithm> // for std::sort() (needed for Tesla machine only)
#include "math.h"
#include "..\Genetic_Algorithm\sort.h"

//Constructor for the output struct
output::output(const cudaConstants* cConstants, const std::string& selectFolder) {
  //Set the base folder to the desired folder
  baseFolder = selectFolder; 

  //Prepare the folders for the files
  prepareFolders(cConstants); 

  //Start the genPreformance file
  initializeGenPerformance(cConstants);

  //Start the simple genPreformance file
  //initializeSimpleGenPerformance(cConstants);
}

//Call the functions which print during a run
void output::printGeneration(const cudaConstants * cConstants, const std::vector<Adult>& allAdults, const std::vector<double>& objectiveAvgValues, const int& generation, const double& new_anneal, int& errorNum, const int& duplicateNum, const int & minSteps, const int & avgSteps, const int & maxSteps, const double& minDist, const double& avgDist, const double& maxDist, const double& avgAge, const int& oldestAge, const double& avgBirthday, const int& oldestBirthday, const float& avgGenTime) {
  
  //Check to see if the best adults should be printed to the terminal on this generation
  if (generation % cConstants->disp_freq == 0) {
    printBestAdults(cConstants, allAdults, generation, errorNum, duplicateNum, oldestBirthday);
  }

  //The next prints should only occur if record mode is on
  if (cConstants->record_mode == true) {

    //Check which sees if this is a write_freq generation
    if (generation % cConstants->write_freq == 0) {

      //Push info to the function which prints the full genPerformance file
      recordGenerationPerformance(cConstants, allAdults, objectiveAvgValues, generation, new_anneal, errorNum, duplicateNum, minSteps, avgSteps, maxSteps, minDist, avgDist, maxDist, avgAge, oldestAge, avgBirthday, oldestBirthday, avgGenTime);
      //Push info to the function which prints the simplified genPerformance file
      //recordGenSimple(cConstants, allAdults, objectiveAvgValues, generation);
    }

    //Check which sees if the full adult list should be printed
    if (generation % cConstants->all_write_freq == 0) {

      //Call the function which prints all adults for this generation into a new file
      recordAllIndividuals("AllAdults", cConstants, allAdults, generation);
    }
  }
}

//Function will handle printing at the end of a run
void output::printFinalGen(const cudaConstants * cConstants, std::vector<Adult>& allAdults, GPUMem & gpuValues, const bool& converged, const int& generation, int& errorNum, const int& duplicateNum, const int& oldestBirthday, const float& avgGenTime) {
  //Call for a print of allIndividuals if in record mode and if it is not a report generation already
  //  Note: second check is simply so redundant files aren't created
  if ((cConstants->record_mode == true) && (generation % cConstants->all_write_freq != 0)) {
    //Print the adults of this final generation
    recordAllIndividuals("AllAdults", cConstants, allAdults, generation);
  }

  //Print the best adults for the last generation to the console
  printBestAdults(cConstants, allAdults, generation, errorNum, duplicateNum, oldestBirthday);

  //Print the result of the run to the combined run result file
  reportRun(cConstants, allAdults, converged, generation, avgGenTime);

  //Check to see if there is a convergence before printing the trajectory
  //if (converged) {
      // Evaluate and print this solution's information to binary files
      trajectoryPrint(generation, cConstants, allAdults[0], gpuValues);
  //}
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Initialize the genPerformance .csv file with header row
void output::initializeGenPerformance(const cudaConstants * cConstants) {
  std::ofstream excelFile;
  int seed = cConstants->time_seed;
  
  // setting the numeric id tag as the randomization seed (when doing runs of various properties, suggested to add other values to differentiate)
  std::string fileId = std::to_string(seed);

  //Create the new genPerformance file
  excelFile.open(outputPath + "genPerformance-" + fileId + ".csv", std::ios_base::app);
  
  excelFile << "gen,";

  //Names for each objective
  for (int i = 0; i < cConstants->missionObjectives.size(); i++) {
    excelFile << "algorithmBest" << cConstants->missionObjectives[i].name << ",";
    excelFile << "algorithmBest" << cConstants->missionObjectives[i].name << "Normalization,";
  }
  for (int i = 0; i < cConstants->missionObjectives.size(); i++) {
    excelFile << "best" << cConstants->missionObjectives[i].name << ",";
    excelFile << "avg" << cConstants->missionObjectives[i].name << ",";
  }

  excelFile << "anneal,avgGenTime,minDistance,avgDistance,maxDistance,minSteps,avgSteps,maxSteps,avgAge,oldestAge,bestAdultAge,avgBirthday,oldestBirthday,bestAdultBirthday,errorNum,duplicateNum,avgParentProgress,progress,parentChildProgressRatio,alpha,beta,zeta,tripTime,";
  
  for (int i = 0; i < GAMMA_ARRAY_SIZE; i++) {
    excelFile << "gamma"; 
    if (i == 0) {
      excelFile << "_a" << i << ",";
    }
    else if (i % 2 == 0) {
      excelFile << "_b" << i/2 << ",";
    }
    else {
      excelFile << "_a" << (i+1)/2 << ",";
    }
  }
  for (int i = 0; i < TAU_ARRAY_SIZE; i++) {
    excelFile << "tau"; 
    if (i == 0) {
      excelFile << "_a" << i << ",";
    }
    else if (i % 2 == 0) {
      excelFile << "_b" << i/2 << ",";
    }
    else {
      excelFile << "_a" << (i+1)/2 << ",";
    }
  }
  for (int i = 0; i < COAST_ARRAY_SIZE; i++) {
    excelFile << "psi"; 
    if (i == 0) {
      excelFile << "_a" << i << ",";
    }
    else if (i % 2 == 0) {
      excelFile << "_b" << i/2 << ",";
    }
    else {
      excelFile << "_a" << (i+1)/2 << ",";
    }
  }

  excelFile << "\n";

  excelFile.close();
}

// Initialize simpleGenPerformance with header rows
void output::initializeSimpleGenPerformance(const cudaConstants* cConstants) {
  
  std::ofstream excelFile;
  int seed = cConstants->time_seed;
  
  // setting the numeric id tag as the randomization seed (when doing runs of various properties, suggested to add other values to differentiate)
  std::string fileId = std::to_string(seed);

  //Create the new genPerformance file
  excelFile.open(outputPath + "simpleGenPerformance-" + fileId + ".csv", std::ios_base::app);

  //Generation and progress
  excelFile << "gen,Rank-Rarity Progress,";

  //Names for each objective

  for (int i = 0; i < cConstants->missionObjectives.size(); i++) {
    excelFile << "rankRarity" << cConstants->missionObjectives[i].name << ",";
  }
  for (int i = 0; i < cConstants->missionObjectives.size(); i++) {
    excelFile << "best" << cConstants->missionObjectives[i].name << ",";
    excelFile << "avg" << cConstants->missionObjectives[i].name << ",";
  }

  //End the line and close the file
  excelFile << "\n";
  excelFile.close();
}

// Checks if needed output files exist, creates them if they don't exist, and copies the run's config files to them
void output::prepareFolders(const cudaConstants* cConstants) {
  //Create output files folder
  if (mkdir(baseFolder.c_str()) == 0) {
    //Report that the output files were created
    std::cout << "\n~~~" << baseFolder << " folder created for the first time~~~\n"; 
  }

  //Version will allow for mutltiple output folders to be created for one seed
  int version = 1; 

  //Flag will determine when a seed folder has been successfully created
  bool seedFolderCreated = false; 

  //While loop will keep trying to make seed folders until a folder has been created or it has tried too many times
  //This will allow for mutliple result folders if using the same seed (to a point)
  while (!seedFolderCreated && version <= cConstants->run_count) {
    
    //Determine the seedFileName string
    //It will use the base folder, then inside the base create a new folder for this run based on the seed and the version
    outputPath = baseFolder + "\\" + std::to_string(static_cast<int>(cConstants->time_seed)) + "_" + std::to_string(version) + "\\";

    //Attempt to open the seed file
    if (mkdir(outputPath.c_str()) == 0) {
      //Seed folder created correctly, setting seed file created to true
      seedFolderCreated = true; 
    }
    else {
      //Seed folder not created correctly, either an issue with the folder, or it already exists
      //Report the potential problem
      std::cout << "\n~~~Potential issue with creating seed folder;\tfolder name tried: " << outputPath << "~~~\n";

      //Assuming it already exists, add one to version to try again with a new version number
      version++;
    }
  }

  //Final check to see if the seed folder was created or not
  //If version is above the run count, it means that the seed folder was not created correctly
  if (version > cConstants->run_count) {
    //Return, since there is no point in trying to copy the config files to the seed folder
    return;
  }

  //Copy the necessary config files into the new folder
  copyFile("..\\Config_Constants\\mission.config", outputPath + "mission.config");
  copyFile("..\\Config_Constants\\genetic.config", outputPath + "genetic.config");
  copyFile("..\\Config_Constants\\" + cConstants->destination, outputPath + cConstants->destination);

  //Return now that the files are prepped
  return; 
}

// Assist function to prepareFolders, which will copy the config files into a seed output folder
void output::copyFile (const std::string& source, const std::string& destination) {
  //Create the file streams for the source and desination
  std::ifstream src(source, std::ios::binary);
  std::ofstream dest(destination, std::ios::binary | std::ios::trunc); 

  //Copy the source to the destination
  dest << src.rdbuf(); 

  //Close the files
  src.close();
  dest.close();
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Take in the current state of the generation and appends to excel file, assumes initializeRecord() had already been called before (no need to output a header row)
void output::recordGenerationPerformance(const cudaConstants * cConstants, std::vector<Adult> adults, const std::vector<double>& objectiveAvgValues, const int& generation, const double& new_anneal, const int& errorNum, const int& duplicateNum, const int & minSteps, const int & avgSteps, const int & maxSteps, const double& minDist, const double& avgDist, const double& maxDist, const double& avgAge, const int& oldestAge, const double& avgBirthday, const int& oldestBirthday, const float& avgGenTime) {
  std::ofstream excelFile;
  int seed = cConstants->time_seed;
  std::string fileId = std::to_string(seed);
  excelFile.open(outputPath + "genPerformance-" + fileId + ".csv", std::ios_base::app);
  excelFile << std::setprecision(20);
  excelFile << std::fixed;
  // Record best individuals best posDiff and speedDiff of this generation
  excelFile << generation << ",";

  //Output the best rank distance adult's parameters
  for (int i = 0; i < cConstants->missionObjectives.size(); i++) {
    if (cConstants->missionObjectives[i].goal < 0) {
      excelFile << adults[0].getParameters(cConstants->missionObjectives[i]) << ",";
    }
    else {
      excelFile << -(adults[0].getParameters(cConstants->missionObjectives[i])) << ",";
    }
    excelFile << adults[0].normalizedObj[i] << ",";
  }
  //Output the objective parameter for the best adult and the average for each objective
  for (int i = 0; i < cConstants->missionObjectives.size(); i++) {
    //Sort the Adults vector by the parameter
    parameterSort(adults, cConstants->missionObjectives[i], adults.size()); 
    
    //Output the 1st adult's value for that parameter, which will be the best of that value in the Adults vector
    //Output the negative if it is a maximization
    if (cConstants->missionObjectives[i].goal < 0) {
      excelFile << adults[0].getParameters(cConstants->missionObjectives[i]) << ",";
    }
    else {
      excelFile << -(adults[0].getParameters(cConstants->missionObjectives[i])) << ",";
    }
    
    //Output the average value for this parameter
    //The 
    excelFile << objectiveAvgValues[i] << ",";
  }

  //Reset the sort to rankRarity
  // std::sort(adults.begin(), adults.end(), rankRaritySort);
  mainSort(adults, cConstants, adults.size());

  //New anneal every gen
  excelFile << new_anneal << ",";

  //Avg generation time
  excelFile << avgGenTime << ",";

  //Distance values
  excelFile << minDist << ",";
  excelFile << avgDist << ",";
  excelFile << maxDist << ",";
  //Step values
  excelFile << minSteps << ",";
  excelFile << avgSteps << ",";
  excelFile << maxSteps << ",";
  //Age values
  excelFile << avgAge << ",";
  excelFile << oldestAge << ",";
  excelFile << generation - adults[0].birthday << ",";
  excelFile << avgBirthday << ",";
  excelFile << oldestBirthday << ",";
  excelFile << adults[0].birthday << ",";
  //Status Counts
  excelFile << errorNum << ",";
  excelFile << duplicateNum << ",";
  //Progress values
  excelFile << adults[0].avgParentProgress << ",";
  excelFile << adults[0].progress << ",";
  excelFile << adults[0].avgParentProgress/adults[0].progress << ",";

  // Record best individual's parameters
  excelFile << adults[0].startParams.alpha << ",";
  excelFile << adults[0].startParams.beta << ",";
  excelFile << adults[0].startParams.zeta << ",";
  excelFile << adults[0].startParams.tripTime << ",";

  for (int i = GAMMA_OFFSET; i < GAMMA_ARRAY_SIZE + GAMMA_OFFSET; i++) {
    excelFile << adults[0].startParams.coeff.gamma[i-GAMMA_OFFSET] << ","; 
  }
  for (int i = TAU_OFFSET; i < TAU_ARRAY_SIZE + TAU_OFFSET; i++) {
    excelFile << adults[0].startParams.coeff.tau[i-TAU_OFFSET] << ","; 
  }
  for (int i = COAST_OFFSET; i < COAST_ARRAY_SIZE + COAST_OFFSET; i++) {
    excelFile << adults[0].startParams.coeff.coast[i-COAST_OFFSET] << ","; 
  } 

  excelFile << "\n"; // End of row
  excelFile.close();
  
}

// Record highlights of the full genPerformance file
void output::recordGenSimple (const cudaConstants* cConstants, std::vector<Adult> adults, const std::vector<double>& objectiveAvgValues, const int& generation) {
  //First open the output file
  std::ofstream excelFile;
  excelFile.open(outputPath + "simpleGenPerformance-" + std::to_string(static_cast<int>(cConstants->time_seed)) + ".csv", std::ios_base::app);

  //Output the generation and the progress of the best rank-distance adult
  excelFile << generation << "," << adults[0].progress << ","; 

  //Output the best rank distance adult's parameters
  for (int i = 0; i < cConstants->missionObjectives.size(); i++) {
    excelFile << adults[0].getParameters(cConstants->missionObjectives[i]) << ",";
  }
  //Output the objective parameter for the best adult and the average for each objective
  for (int i = 0; i < cConstants->missionObjectives.size(); i++) {
    //Sort the Adults vector by the parameter
    parameterSort(adults, cConstants->missionObjectives[i], adults.size()); 
    
    //Output the 1st adult's value for that parameter, which will be the best of that value in the Adults vector
    excelFile << adults[0].getParameters(cConstants->missionObjectives[i]) << ",";
    //Output the average value for this parameter
    excelFile << objectiveAvgValues[i] << ",";
  }

  //Re-sort the adults in rank-rarity order
  // std::sort(adults.begin(), adults.end(), rankRaritySort);
  mainSort(adults, cConstants, adults.size());

  //End the line for this generation and close the file
  excelFile << "\n"; 
  excelFile.close(); 
}

// Takes in an adult vector and records the parameter info on all individuals within the vector
void output::recordAllIndividuals(std::string name, const cudaConstants * cConstants, const std::vector<Adult>& adults, const int& generation) {
  
  std::ofstream outputFile;
  outputFile.open(outputPath + std::to_string(static_cast<int>(cConstants->time_seed)) + "-" + name +"-gen#" + std::to_string(generation) + ".csv");

  // Setup the header row
  outputFile << "position,";

  for (int i = 0; i < cConstants->missionObjectives.size(); i++) {
    outputFile << cConstants->missionObjectives[i].name << ",";
    outputFile << cConstants->missionObjectives[i].name << "Normalization,";
  }
  
  outputFile << "rank,rarity,distance,numSteps,age,birthday,avgParentProgress,progress,parentChildProgressRatio,";

  for (int i = 0; i < GAMMA_ARRAY_SIZE; i++) {
    outputFile << "gamma" << i << ",";
  }
  for (int i = 0; i < TAU_ARRAY_SIZE; i++) {
    outputFile << "tau" << i << ",";
  }
  outputFile << "alpha,beta,zeta,tripTime,";
  for (int i = 0; i < COAST_ARRAY_SIZE; i++) {
    outputFile << "coast" << i << ",";
  }

  outputFile << '\n';

  outputFile << std::setprecision(20);
  outputFile << std::fixed;

  // Record all individuals in the adults vector
  for (int i = 0; i < adults.size(); i++) {
    outputFile << i << ",";

    for (int j = 0; j < cConstants->missionObjectives.size(); j++)
    {
      outputFile << adults[i].getParameters(cConstants->missionObjectives[j]) << ",";
      outputFile << adults[i].normalizedObj[j] << ",";
    }
    
    outputFile << adults[i].rank << ",";
    outputFile << adults[i].rarity << ",";
    outputFile << adults[i].distance << ",";
    outputFile << adults[i].stepCount << ",";
    outputFile << generation-adults[i].birthday << ",";
    outputFile << adults[i].birthday << ",";
    outputFile << adults[i].avgParentProgress << ",";
    outputFile << adults[i].progress << ",";
    outputFile << adults[i].progress/adults[i].avgParentProgress << ",";

    for (int j = 0; j < GAMMA_ARRAY_SIZE; j++) {
      outputFile << adults[i].startParams.coeff.gamma[j] << ",";
    }
    for (int j = 0; j < TAU_ARRAY_SIZE; j++) {
      outputFile << adults[i].startParams.coeff.tau[j] << ",";
    }
    outputFile << adults[i].startParams.alpha << ",";
    outputFile << adults[i].startParams.beta << ",";
    outputFile << adults[i].startParams.zeta << ",";
    outputFile << adults[i].startParams.tripTime << ",";
    for (int j = 0; j < COAST_ARRAY_SIZE; j++) {
      outputFile << adults[i].startParams.coeff.coast[j] << ",";
    }
    outputFile << "\n";
  }
  outputFile.close();
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Function which prints the general info of a run when its over
void output::reportRun(const cudaConstants* cConstants, const std::vector<Adult>& adults, const bool& converged, const int& generation, const float& avgGenTime) {
  //Open the runReport file
  std::ofstream output("..\\Output_Files\\runReports.csv", std::ios::app);

  //Report the seed
  output << "Seed:," << static_cast<int>(cConstants->time_seed) << ",,";
  //Report the final generion
  output << "Final Generation:," << generation << ",,";
  //Report the average generation time
  output << "Average Time/Generation:," << avgGenTime << ",,";
  //Report if it converged
  output << "Converged:,";
  if (converged) {
    output << "yes\n";
  }
  else {
    output << "no\n";
  }

  //Create a key for the parameters, progress, and age 
  //Output basic adult stats
  output << "R-D Position,Progress,Age";
  //First output parameters
  for (int i = 0; i < cConstants->missionObjectives.size(); i++) {
    //Output the objective name
    output << "," << cConstants->missionObjectives[i].name;
  }
  //End line to move on to reporting adults
  output << "\n";

  //Report the top 3 rank-distance adults
  for (int i = 0; i < 3; i++){
    //Report their position, progress, and age
    output << i << "," << adults[i].progress << "," << generation - adults[i].birthday; 

    //For each objective, print the adult's parameter
    for (int j = 0; j < cConstants->missionObjectives.size(); j++) {
      output << "," << adults[i].getParameters(cConstants->missionObjectives[j]);
    }

    //End line to move on to the next adult
    output << "\n";
  }
  //End line before closing to make space for the next run
  output << "\n";
  output.close();
}

// Main Output, final results of genetic algorithm
void output::trajectoryPrint(int generation, const cudaConstants* cConstants, const Adult& best, GPUMem & gpuValues) {
  /*set the target and inital conditions for the earth and spacecraft:
  constructor takes in radial position(au), angluar position(rad), axial position(au),
  radial velocity(au/s), tangential velocity(au/s), axial velocity(au/s)*/

  // setting landing conditions of the target (Sept 30, 2022)
  elements<double> target = elements<double>(cConstants->r_fin_target, cConstants->theta_fin_target, cConstants->z_fin_target, cConstants->vr_fin_target, cConstants->vtheta_fin_target, cConstants->vz_fin_target);

// setting the final position of Mars based on the landing date
  elements<double> mars = elements<double>(cConstants->r_fin_mars, cConstants->theta_fin_mars, cConstants->z_fin_mars, cConstants->vr_fin_mars, cConstants->vtheta_fin_mars, cConstants->vz_fin_mars);

  // setting initial conditions of earth based off of the impact date (Sept 30, 2022) minus the trip time (optimized).
  elements<double> earth =  launchCon->getCondition(best.startParams.tripTime);
  
  // setting initial conditions of the spacecraft
  elements<double> spaceCraft = elements<double>(earth.r+ESOI*cos(best.startParams.alpha),
                                                 earth.theta + asin(sin(M_PI-best.startParams.alpha)*ESOI/earth.r),
                                                 earth.z,
                                                 earth.vr + cos(best.startParams.zeta)*sin(best.startParams.beta)*cConstants->v_escape,
                                                 earth.vtheta + cos(best.startParams.zeta)*cos(best.startParams.beta)*cConstants->v_escape,
                                                 earth.vz + sin(best.startParams.zeta)*cConstants->v_escape);

  //Calculate the maximum possible number of steps taken to make enough memory is allocated
  int numSteps = (cConstants->maxSimNum * (cConstants->max_numsteps + 1))+1;

  // Assigning wetMass
  double wetMass = cConstants->wet_mass;

  // Initialize memory for the solution vector of the dependant solution
  elements<double>* yp = new elements<double>[numSteps];
  elements<double>* marsIndex = new elements<double>[numSteps];
  
  double *times, *gamma, *tau, *accel_output, *fuelSpent, *work, *dE, *Etot_avg;
  times = new double[numSteps]; // Initialize memory for time array
  gamma = new double[numSteps]; // Initialize memory for gamma array
  tau = new double[numSteps]; // Initialize memory for tau array
  accel_output = new double[numSteps]; // Initialize memory for acceleration array
  fuelSpent = new double[numSteps];  // Initialize memory for fuelSpent array
  work = new double[numSteps];  // Initialize memory for work array
  dE = new double[numSteps];  // Initialize memory for delta-E array
  Etot_avg = new double[numSteps];  // Initialize memory for average mechanical energy array

  //store the number of threads
  const int numThreads = gpuValues.numThreads;
  const int blockThreads = cConstants->thread_block_size;

  //Create new dummy child from the best Adult
  rkParameters<double> bestParams = best.startParams;
  Child *convergedChild = new Child(bestParams, cConstants, best.birthday, best.avgParentProgress);

  //Will simulate the best individual until they have completed their simulation
  while ((*convergedChild).simStatus != COMPLETED_SIM) {

    //Create CUDA event
    cudaEvent_t kernelStart, kernelEnd;
    cudaEventCreate(&kernelStart);
    cudaEventCreate(&kernelEnd);

    // copy values of parameters passed from host onto device
    cudaMemcpy(gpuValues.devGeneration, convergedChild, sizeof(Child), cudaMemcpyHostToDevice);

    //Run CUDA simulation
    cudaEventRecord(kernelStart);
    rk4CUDASim<<<(numThreads+blockThreads-1)/blockThreads,blockThreads>>>(gpuValues.devGeneration, gpuValues.devAbsTol, 1, gpuValues.devCConstant, gpuValues.devMarsLaunchCon, gpuValues.devTime_steps, gpuValues.devY_steps, gpuValues.devGamma_steps, gpuValues.devTau_steps, gpuValues.devAccel_steps, gpuValues.devFuel_steps);
    cudaEventRecord(kernelEnd);

    // copy the result of the kernel onto the host
    cudaMemcpy(convergedChild, gpuValues.devGeneration, sizeof(Child), cudaMemcpyDeviceToHost);

    //wait for CUDA sim to finish
    float kernelT;   
    cudaEventSynchronize(kernelEnd);
    cudaEventElapsedTime(&kernelT, kernelStart, kernelEnd);
  }

  //Copy the step-by-step values off of the GPU
  cudaMemcpy(times, gpuValues.devTime_steps, numSteps * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(yp, gpuValues.devY_steps, numSteps * sizeof(elements<double>), cudaMemcpyDeviceToHost);
  cudaMemcpy(gamma, gpuValues.devGamma_steps, numSteps * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(tau, gpuValues.devTau_steps, numSteps * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(accel_output, gpuValues.devAccel_steps, numSteps * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(fuelSpent, gpuValues.devFuel_steps, numSteps * sizeof(double), cudaMemcpyDeviceToHost);

  // calculate the error in conservation of mechanical energy due to the thruster
  errorCheck(times, yp, gamma, tau, static_cast<int>((*convergedChild).stepCount), accel_output, fuelSpent, wetMass, work, dE, Etot_avg, cConstants, marsIndex);

  //Get the seed for outputs
  int seed = cConstants->time_seed;

  // Test outputs
  std::cout << "\nComparison\n";
  for (int i = 0; i < cConstants->missionObjectives.size(); i++) {
      std::cout << "Main Run " << cConstants->missionObjectives[i].name << ": ";
      std::cout << best.getParameters(cConstants->missionObjectives[i]) << "\n";
  }
  std::cout << "Main Run steps: " << best.stepCount << "\n\n";

  //Get equivalent outputs for convergedChild
  (*convergedChild).getPosDiff(cConstants);
  (*convergedChild).getSpeedDiff(cConstants);
  (*convergedChild).getOrbitPosDiff(cConstants);
  (*convergedChild).getOrbitSpeedDiff(cConstants);

  //Outputting equivalent outputs
  for (int i = 0; i < cConstants->missionObjectives.size(); i++) {
      std::cout << "Bin Run " << cConstants->missionObjectives[i].name << ": ";
      std::cout << (*convergedChild).getParameters(cConstants->missionObjectives[i]) << "\n";
  }
  std::cout << "Bin Run steps: " << (*convergedChild).stepCount << "\n\n";

  // OUTDATED SINCE 2023
  // This function is used to compare the final best thread with other runs
  // append this thread's info to a csv file
  //if (cConstants->record_mode == true) {
    // Record initial and final fuel masses along with tripTime and relative velocity at impact
    //recordFuelOutput(cConstants, x, fuelSpent[lastStepInt], best, seed);
    //progressiveAnalysis(generation, lastStepInt, x, yOut, cConstants);
  //}

  // binary outputs
  std::ofstream output;

  //step-by-step bin file
  output.open(outputPath + "orbitalMotion-"+std::to_string(seed)+".bin", std::ios::binary);

  // Output this thread's data at each time step
  for(int i = 0; i <= (*convergedChild).stepCount; i++) {
    output.write((char*)&yp[i], sizeof (elements<double>));
    output.write((char*)&times[i], sizeof (double));
    output.write((char*)&gamma[i], sizeof (double));
    output.write((char*)&tau[i], sizeof (double));
    output.write((char*)&accel_output[i], sizeof (double));
    output.write((char*)&fuelSpent[i], sizeof (double));
    output.write((char*)&work[i], sizeof(double));
    output.write((char*)&dE[i], sizeof(double));
    output.write((char*)&Etot_avg[i], sizeof(double));
  }
  output.close();


  //Setup conditions file
  output.open(outputPath + "finalOptimization-" + std::to_string(seed)+".bin", std::ios::binary);

  double gsize = GAMMA_ARRAY_SIZE, tsize = TAU_ARRAY_SIZE, csize = COAST_ARRAY_SIZE;

  // Impact conditions
  output.write((char*)&cConstants->r_fin_target, sizeof(double));
  output.write((char*)&cConstants->theta_fin_target, sizeof(double));
  output.write((char*)&cConstants->z_fin_target, sizeof(double));
  output.write((char*)&cConstants->vr_fin_target, sizeof(double));
  output.write((char*)&cConstants->vtheta_fin_target, sizeof(double));
  output.write((char*)&cConstants->vz_fin_target, sizeof(double));

  output.write((char*)&cConstants->r_fin_mars, sizeof(double));
  output.write((char*)&cConstants->theta_fin_mars, sizeof(double));
  output.write((char*)&cConstants->z_fin_mars, sizeof(double));
  output.write((char*)&cConstants->vr_fin_mars, sizeof(double));
  output.write((char*)&cConstants->vtheta_fin_mars, sizeof(double));
  output.write((char*)&cConstants->vz_fin_mars, sizeof(double));

  output.write((char*)&cConstants->r_fin_earth, sizeof(double));
  output.write((char*)&cConstants->theta_fin_earth, sizeof(double));
  output.write((char*)&cConstants->z_fin_earth, sizeof(double));
  output.write((char*)&cConstants->vr_fin_earth, sizeof(double));
  output.write((char*)&cConstants->vtheta_fin_earth, sizeof(double));
  output.write((char*)&cConstants->vz_fin_earth, sizeof(double));

  // Launch conditions
  output.write((char*)&earth.r, sizeof(double));
  output.write((char*)&earth.theta, sizeof(double));
  output.write((char*)&earth.z, sizeof(double));
  
  // Thruster info
  output.write((char*)&cConstants->fuel_mass, sizeof(double));
  output.write((char*)&cConstants->coast_threshold, sizeof(double));
  output.write((char*)&gsize, sizeof(double));
  output.write((char*)&tsize, sizeof(double));
  output.write((char*)&csize, sizeof(double));

  // Optimized starting parameters
  for (int i = 0; i < GAMMA_ARRAY_SIZE; i++) {
    output.write((char*)&best.startParams.coeff.gamma[i], sizeof (double));
  }
  for (int i = 0; i < TAU_ARRAY_SIZE; i++) {
    output.write((char*)&best.startParams.coeff.tau[i], sizeof (double)); 
  }
  output.write((char*)&best.startParams.alpha, sizeof (double));
  output.write((char*)&best.startParams.beta, sizeof (double));
  output.write((char*)&best.startParams.zeta, sizeof (double));
  output.write((char*)&best.startParams.tripTime, sizeof (double));
  for (int i = 0; i < COAST_ARRAY_SIZE; i++) {
    output.write((char*)&best.startParams.coeff.coast[i], sizeof (double));
  }

  //Matlab can only take in one data type when pulling in the data, so beacuse everything else is a double, step counts need to be too
  double lastStepDouble = static_cast<double>((*convergedChild).stepCount);
  
  // Number of steps taken in final RK calculation
  output.write((char*)&lastStepDouble, sizeof (double));

  output.close();
  
  // cleaning up dynamic memory
  delete [] yp;
  delete [] times;
  delete [] gamma;
  delete [] tau;
  delete [] accel_output;
  delete [] fuelSpent;
  delete [] work;
  delete [] dE;
  delete [] Etot_avg;
}

// Records error in energy conservation due to thrust calculations
void output::errorCheck(double *time, elements<double> *yp,  double *gamma,  double *tau, int & lastStep, double *accel, double *fuelSpent, const double & wetMass, double *work, double *dE, double *Etot_avg, const cudaConstants* config, elements<double> *mars) {
  // Initializing storage at each time step
  double *mass = new double[lastStep];
  double *Etot = new double[lastStep];
  double marsCraftDist;

  //Iterate over all time steps
  for (int i = 0; i < lastStep; i++) {
    // total mass at each time step
    mass[i] = wetMass - fuelSpent[i];
    // total mechanical energy (K + U) at each time step
    marsCraftDist = sqrt(pow(mars[i].r, 2) + pow(yp[i].r, 2) + pow(yp[i].z - mars[i].z, 2) - (2*yp[i].r*mars[i].r*cos(mars[i].theta-yp[i].theta)));

    Etot[i] = mass[i] * ((pow(yp[i].vr,2) + pow(yp[i].vtheta,2) + pow(yp[i].vz,2))/ 2 - (constG * massSun / yp[i].r) - (constG * massMars/ marsCraftDist )) / pow(AU,2);
    if (i) {
      // Ignore first time step, only calculate as a change from step to step
      // W = F.dL (work done between time steps)
      work[i] = (mass[i] + mass[i-1])/2 * (accel[i] + accel[i-1])/2 * ((sin((gamma[i] + gamma[i-1])/2)*cos((tau[i] + tau[i-1])/2)*(yp[i].r - yp[i-1].r)) + (cos((gamma[i] + gamma[i-1])/2)*cos((tau[i] + tau[i-1])/2)*(yp[i].r + yp[i-1].r)/2*(yp[i].theta - yp[i-1].theta)) + (sin((tau[i] + tau[i-1])/2)*(yp[i].z - yp[i-1].z))) / pow(AU,2);
      // change in mechanical energy between time steps
      dE[i] = Etot[i] - Etot[i-1];
      // average mechanical energy between time steps
      Etot_avg[i] = (Etot[i] + Etot[i-1])/2;
    }
  }
  work[0] = dE[0] = 0;
  Etot_avg[0] = Etot[0];

  // Output data is passed by reference into trajectoryPrint and integrated into orbitalMotion-seed.bin

  // cleaning up dynamic memory
  delete [] mass;
  delete [] Etot;
}

void output::recordEarthData(const cudaConstants * cConstants, const int & generation) {
  double timeStamp = cConstants->triptime_min; // Start the timeStamp at the triptime min (closest to impact)
  // seed to hold time_seed value to identify the file
  int seed = cConstants->time_seed;
  // Open file
  std::ofstream planetValues;
  planetValues.open(outputPath + "earthCheckValues-"+ std::to_string(seed)+ "-gen#" + std::to_string(generation) + ".csv");
  // Set header row for the table to record values, with timeStamp
  planetValues << "TimeStamp, Radius, Theta, Z, vRadius, vTheta, vZ\n";
  // Fill in rows for elements of launchCon across the time range going backwards until triptime max
  while (timeStamp < cConstants->triptime_max) {
      planetValues << timeStamp << "," << launchCon->getCondition(timeStamp);
      timeStamp += cConstants->timeRes; // Increment to next hour as timeRes is assumed to be set to every hour in this output
  }
  
  // Done recording earth calculations, close file
  planetValues.close();
}

// Stores information of launchCon of timeRes (was timeRes*24) resolution in EarthCheckValues-[time_seed].csv
void output::recordMarsData(const cudaConstants * cConstants, const int & generation) {
  double timeStamp = 0.0; // Start the timeStamp at the triptime min (closest to impact)
  // seed to hold time_seed value to identify the file
  int seed = cConstants->time_seed;
  // Open file
  std::ofstream planetValues;
  planetValues.open(outputPath + "marsCheckValues-"+ std::to_string(seed)+ "-gen#" + std::to_string(generation) + ".csv");
  // Set header row for the table to record values, with timeStamp
  planetValues << "TimeStamp, Radius, Theta, Z, vRadius, vTheta, vZ\n";
  // Fill in rows for elements of launchCon across the time range going backwards until triptime max
  while (timeStamp < cConstants->triptime_max) {
      planetValues << timeStamp << "," << marsLaunchCon->getCondition(timeStamp);
      timeStamp += cConstants->timeRes; // Increment to next hour as timeRes is assumed to be set to every hour in this output
  }
  
  // Done recording earth calculations, close file
  planetValues.close();
}

//Method used to print a run's used reference points
void output::recordReferencePoints(const cudaConstants * cConstants, const ReferencePoints & refPoints) {

  // seed to hold time_seed value to identify the file
  int seed = cConstants->time_seed;

  // Open file
  std::ofstream output;
  output.open(outputPath + "referencePoints-"+ std::to_string(seed) + ".csv");

  //Print the number of objectives as a header
  for (int i = 0; i < refPoints.points[0].size(); i++) {
    output << "Objective" << i << ",";  
  }

  //Print all of the points
  for (int i = 0; i < refPoints.points.size(); i++) {
    //New line for new point
    output << "\n";

    for (int j = 0; j < refPoints.points[i].size(); j++) {

      //Print one part of the point
      output << refPoints.points[i][j] << ",";
    }
  }
  
  //Close the file
  output.close();
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//General Functions

//Function which will print the best rank-distance adult and the best adult for each objective to the terminal
void printBestAdults(const cudaConstants* cConstants, std::vector<Adult> adults, const int& generation, int& numErrors, const int& numDuplicates, const int& oldestBirthday) {
  // Prints the best individual's posDiff / speedDiff

  //Print the generation
  std::cout << "\n\nGeneration " << generation << " data:\n";

  //Print the best adult for each objective
  for (int i = 0; i < cConstants->missionObjectives.size(); i++) {
      //Sort the adults array by the correct parameter & order
      parameterSort(adults, cConstants->missionObjectives[i], adults.size());

      //Print the name of the objective
      std::cout << "\nBest " << cConstants->missionObjectives[i].name << " Individual:";
      //Display the info to the terminal
      terminalDisplay(adults[0], cConstants->missionObjectives, cConstants);
  }

  //display to the terminal the best individual based on rankRarity
  std::cout << "\nBest Overall Individual:";

  // std::sort(adults.begin(), adults.end(), rankRaritySort);
  mainSort(adults, cConstants, adults.size());

  terminalDisplay(adults[0], cConstants->missionObjectives, cConstants);

  //Display number of errors
  std::cout << "\n# of errors this generation: " << numErrors << "\n";

  //Display number of duplicates
  std::cout << "\n# of duplicates this generation: " << numDuplicates << "\n";
  
  //display the oldest individual
  std::cout << "\nOldest age adult: " << generation - oldestBirthday << "\n";

  //Display the steps taken by the best individual
  std::cout << "\nBest step count: " << adults[0].stepCount << "\n";

  //Display the progress of the best rank distance individual 
  std::cout << "\nBest rank-rarity adult progress: " << adults[0].progress << "\n\n";
}

// Utility function to display the currently best individual onto the terminal while the algorithm is still running
void terminalDisplay(const Adult& individual, const std::vector<objective> objectives, const cudaConstants* cConstants) {
    
  //Print the parameters for each of the objectives for the passed in individual
  for (int i = 0; i < objectives.size(); i++) {
    //Print the name and value of the data of the objective
    if (objectives[i].goal < 0){
      std::cout << "\n\t" << objectives[i].name << ": " << individual.getParameters(objectives[i]);
    }
    else {
      std::cout << "\n\t" << objectives[i].name << ": " << -(individual.getParameters(objectives[i]));
    }

    //Print the normalization of the objective if the run is using rank-rarity
    if (cConstants->algorithm == RANK_RARITY) {
      std::cout << "\n\t" << objectives[i].name << " normalization: " << individual.normalizedObj[i];
    }
    
  }

  //Print the distance of the adult if the run is using rank-distance
  if (cConstants->algorithm == RANK_DISTANCE) {
     std::cout << "\n\tDistance: " << individual.distance;
  }

  std::cout << std::endl;
}