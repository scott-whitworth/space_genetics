#include <string>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <direct.h>
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
  initializeSimpleGenPerformance(cConstants);
}

//Call the functions which print during a run
void output::printGeneration(const cudaConstants * cConstants, const std::vector<Adult>& allAdults, const std::vector<double>& objectiveAvgValues, const int& generation, const double& new_anneal, int& errorNum, const int& duplicateNum, const double& minDist, const double& avgDist, const double& maxDist, const double& avgAge, const int& oldestAge, const double& avgBirthday, const int& oldestBirthday) {
  
  //Check to see if the best adults should be printed to the terminal on this generation
  if (generation % cConstants->disp_freq == 0) {
    printBestAdults(cConstants, allAdults, generation, errorNum, duplicateNum, oldestBirthday);
  }

  //The next prints should only occur if record mode is on
  if (cConstants->record_mode == true) {

    //Check which sees if this is a write_freq generation
    if (generation % cConstants->write_freq == 0) {

      //Push info to the function which prints the full genPerformance file
      recordGenerationPerformance(cConstants, allAdults, objectiveAvgValues, generation, new_anneal, errorNum, duplicateNum, minDist, avgDist, maxDist, avgAge, oldestAge, avgBirthday, oldestBirthday);
      //Push info to the function which prints the simplified genPerformance file
      recordGenSimple(cConstants, allAdults, objectiveAvgValues, generation);
    }

    //Check which sees if the full adult list should be printed
    if (generation % cConstants->all_write_freq == 0) {

      //Call the function which prints all adults for this generation into a new file
      recordAllIndividuals("AllAdults", cConstants, allAdults, generation);
    }
  }
}

//Function will handle printing at the end of a run
void output::printFinalGen(const cudaConstants * cConstants, std::vector<Adult>& allAdults, const bool& converged, const int& generation, int& errorNum, const int& duplicateNum, const int& oldestBirthday) {
  //Call for a print of allIndividuals if in record mode and if it is not a report generation already
  //  Note: second check is simply so redundant files aren't created
  if ((cConstants->record_mode == true) && (generation % cConstants->all_write_freq != 0)) {
    //Print the adults of this final generation
    recordAllIndividuals("AllAdults", cConstants, allAdults, generation);
  }

  //Print the best adults for the last generation to the console
  printBestAdults(cConstants, allAdults, generation, errorNum, duplicateNum, oldestBirthday);

  //Print the result of the run to the combined run result file
  reportRun(cConstants, allAdults, converged, generation);

  //Check to see if there is a convergence before printing the trajectory
  //if (converged) {
    //Create the trajectory bin file
    //std::sort(allAdults.begin(), allAdults.end(), bestProgress); //remove if not needed, puts the best progress individual first
    finalRecord(cConstants, allAdults[0], generation);
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
    excelFile << "rankDistance" << cConstants->missionObjectives[i].name << ",";
  }
  for (int i = 0; i < cConstants->missionObjectives.size(); i++) {
    excelFile << "best" << cConstants->missionObjectives[i].name << ",";
    excelFile << "avg" << cConstants->missionObjectives[i].name << ",";
  }

  excelFile << "alpha,beta,zeta,tripTime,";
  
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

  excelFile << "anneal,minDistance,avgDistance,maxDistance,avgAge,oldestAge,bestAdultAge,avgBirthday,oldestBirthday,bestAdultBirthday,errorNum,duplicateNum,avgParentProgress,progress,parentChildProgressRatio\n";
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
  excelFile << "gen,Rank-Distance Progress,";

  //Names for each objective

  for (int i = 0; i < cConstants->missionObjectives.size(); i++) {
    excelFile << "rankDistance" << cConstants->missionObjectives[i].name << ",";
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
void output::recordGenerationPerformance(const cudaConstants * cConstants, std::vector<Adult> adults, const std::vector<double>& objectiveAvgValues, const int& generation, const double& new_anneal, const int& errorNum, const int& duplicateNum, const double& minDist, const double& avgDist, const double& maxDist, const double& avgAge, const int& oldestAge, const double& avgBirthday, const int& oldestBirthday) {
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

  //New anneal every gen
  excelFile << new_anneal << ",";

  //Distance values
  excelFile << minDist << ",";
  excelFile << avgDist << ",";
  excelFile << maxDist << ",";
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

  //End the line for this generation and close the file
  excelFile << "\n"; 
  excelFile.close(); 
}

// Takes in an adult vector and records the parameter info on all individuals within the vector
void output::recordAllIndividuals(std::string name, const cudaConstants * cConstants, const std::vector<Adult>& adults, const int& generation) {
  
  std::ofstream outputFile;
  outputFile.open(outputPath + std::to_string(static_cast<int>(cConstants->time_seed)) + "-" + name +"-gen#" + std::to_string(generation) + ".csv");

  // Setup the header row
  outputFile << "position,alpha,beta,zeta,tripTime,";
  for (int i = 0; i < GAMMA_ARRAY_SIZE; i++) {
    outputFile << "gamma" << i << ",";
  }
  for (int i = 0; i < TAU_ARRAY_SIZE; i++) {
    outputFile << "tau" << i << ",";
  }
  for (int i = 0; i < COAST_ARRAY_SIZE; i++) {
    outputFile << "coast" << i << ",";
  }
  for (int i = 0; i < cConstants->missionObjectives.size(); i++) {
    outputFile << cConstants->missionObjectives[i].name << ",";
  }
  
  outputFile << "age,birthday,rank,distance,avgParentProgress,progress,parentChildProgressRatio";
  outputFile << '\n';

  outputFile << std::setprecision(20);
  outputFile << std::fixed;

  // Record all individuals in the adults vector
  for (int i = 0; i < adults.size(); i++) {
    outputFile << i << ",";
    outputFile << adults[i].startParams.alpha << ",";
    outputFile << adults[i].startParams.beta << ",";
    outputFile << adults[i].startParams.zeta << ",";
    outputFile << adults[i].startParams.tripTime << ",";
    

    for (int j = 0; j < GAMMA_ARRAY_SIZE; j++) {
      outputFile << adults[i].startParams.coeff.gamma[j] << ",";
    }
    for (int j = 0; j < TAU_ARRAY_SIZE; j++) {
      outputFile << adults[i].startParams.coeff.tau[j] << ",";
    }
    for (int j = 0; j < COAST_ARRAY_SIZE; j++) {
      outputFile << adults[i].startParams.coeff.coast[j] << ",";
    }
    for (int j = 0; j < cConstants->missionObjectives.size(); j++)
    {
      outputFile << adults[i].getParameters(cConstants->missionObjectives[j]) << ",";
    }
    
    outputFile << generation-adults[i].birthday << ",";
    outputFile << adults[i].birthday << ",";
    outputFile << adults[i].rank << ",";
    outputFile << adults[i].distance << ",";
    outputFile << adults[i].avgParentProgress << ",";
    outputFile << adults[i].progress << ",";
    outputFile << adults[i].avgParentProgress/adults[i].progress << ",";
    outputFile << "\n";
  }
  outputFile.close();
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Method for doing recording information at the end of the optimization process
void output::finalRecord(const cudaConstants* cConstants, const Adult& bestAdult, const int& generation) {
  // To store parameter values and pass onto writeTrajectoryToFile
  double *start = new double[OPTIM_VARS];

  // Output the final best individual
  for (int j = 0; j < bestAdult.startParams.coeff.gammaSize; j++) {
      start[GAMMA_OFFSET + j] = bestAdult.startParams.coeff.gamma[j];
  }
  for (int j = 0; j < bestAdult.startParams.coeff.tauSize; j++) {
      start[TAU_OFFSET + j] = bestAdult.startParams.coeff.tau[j];
  }
  for (int j = 0; j < bestAdult.startParams.coeff.coastSize; j++) {
      start[COAST_OFFSET + j] = bestAdult.startParams.coeff.coast[j];
  }

  start[TRIPTIME_OFFSET] = bestAdult.startParams.tripTime;
  start[ALPHA_OFFSET] = bestAdult.startParams.alpha;
  start[BETA_OFFSET] = bestAdult.startParams.beta;
  start[ZETA_OFFSET] = bestAdult.startParams.zeta;

  // Test outputs
  std::cout << "Comparison\n";
  for (int i = 0; i < cConstants->missionObjectives.size(); i++) {
      std::cout << "CUDA " << cConstants->missionObjectives[i].name << ": ";
      std::cout << bestAdult.getParameters(cConstants->missionObjectives[i]) << "\n";
  }

  // Evaluate and print this solution's information to binary files
  trajectoryPrint(start, generation, cConstants, bestAdult);

  // cleaning up dynamic memory
  delete [] start;
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Function which prints the general info of a run when its over
void output::reportRun(const cudaConstants* cConstants, const std::vector<Adult>& adults, const bool& converged, const int& generation) {
  //Open the runReport file
  std::ofstream output("..\\Output_Files\\runReports.csv", std::ios::app);

  //Report the seed
  output << "Seed:," << static_cast<int>(cConstants->time_seed) << ",,";
  //Report the final generion
  output << "Final Generation:," << generation << ",,";
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
void output::trajectoryPrint( double x[], int generation, const cudaConstants* cConstants, const Adult& best) {
  /*set the target and inital conditions for the earth and spacecraft:
  constructor takes in radial position(au), angluar position(rad), axial position(au),
  radial velocity(au/s), tangential velocity(au/s), axial velocity(au/s)*/

  // setting landing conditions of the target (Sept 30, 2022)
  elements<double> target = elements<double>(cConstants->r_fin_target, cConstants->theta_fin_target, cConstants->z_fin_target, cConstants->vr_fin_target, cConstants->vtheta_fin_target, cConstants->vz_fin_target);

// setting the final position of Mars based on the landing date
  elements<double> mars = elements<double>(cConstants->r_fin_mars, cConstants->theta_fin_mars, cConstants->z_fin_mars, cConstants->vr_fin_mars, cConstants->vtheta_fin_mars, cConstants->vz_fin_mars);

  // setting initial conditions of earth based off of the impact date (Sept 30, 2022) minus the trip time (optimized).
  elements<double> earth =  launchCon->getCondition(x[TRIPTIME_OFFSET]);
  
  // setting initial conditions of the spacecraft
  elements<double> spaceCraft = elements<double>(earth.r+ESOI*cos(x[ALPHA_OFFSET]),
                                                 earth.theta + asin(sin(M_PI-x[ALPHA_OFFSET])*ESOI/earth.r),
                                                 earth.z,
                                                 earth.vr + cos(x[ZETA_OFFSET])*sin(x[BETA_OFFSET])*cConstants->v_escape,
                                                 earth.vtheta + cos(x[ZETA_OFFSET])*cos(x[BETA_OFFSET])*cConstants->v_escape,
                                                 earth.vz + sin(x[ZETA_OFFSET])*cConstants->v_escape);

  // setting time parameters
  double timeInitial=0; 
  double timeFinal=cConstants->triptime_min; // The shortest possible triptime to make the initial deltaT a small but reasonable size
  double deltaT; // time step
  //int numSteps = cConstants->max_numsteps*7; // initial guess for the number of time steps, guess for the memory allocated 
  int numSteps = best.stepCount+5; // initial guess for the number of time steps, guess for the memory allocated 
  deltaT = (timeFinal-timeInitial) / cConstants->max_numsteps; // initial guess for time step, small is preferable

  // setup of thrust angle calculations based off of optimized coefficients
  coefficients<double> coeff;
  initCoefficient(x,coeff, cConstants);
  // Assigning coast threshold (now done in coefficients because is a constant)

  // Assigning wetMass
  double wetMass = cConstants->wet_mass;
  // setting Runge-Kutta tolerance
  double absTol = cConstants->rk_tol;
  // set optmization minimum
  // double Fmin = cConstants->f_min;

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

  // used to track the cost function throughout a run via output and outputs to a binary
  int lastStepInt;

  // integrate the trajectory of the input starting conditions
  rk4sys(timeInitial, x[TRIPTIME_OFFSET] , times, spaceCraft, deltaT, yp, absTol, coeff, gamma, tau, lastStepInt, accel_output, fuelSpent, wetMass, cConstants, marsIndex);

  // store the number of steps as a double for binary output
  double lastStep = lastStepInt;

  // gets the final y values of the spacecrafts for the cost function.
  elements<double> yOut = yp[lastStepInt];

  // calculate the error in conservation of mechanical energy due to the thruster
  errorCheck(times, yp, gamma, tau, lastStepInt, accel_output, fuelSpent, wetMass, work, dE, Etot_avg, cConstants, marsIndex);


  //Get the seed for outputs
  int seed = cConstants->time_seed;

  // This function is used to compare the final best thread with other runs
  // append this thread's info to a csv file
  //if (cConstants->record_mode == true) {
    // Record initial and final fuel masses along with tripTime and relative velocity at impact
    //recordFuelOutput(cConstants, x, fuelSpent[lastStepInt], best, seed);
    //progressiveAnalysis(generation, lastStepInt, x, yOut, cConstants);
  //}

  // binary outputs
  std::ofstream output;
  
  output.open(outputPath + "orbitalMotion-"+std::to_string(seed)+".bin", std::ios::binary);
  // output.open("orbitalMotion-"+std::to_string(static_cast<int>(seed))+"-"+std::to_string(threadRank)+".bin", std::ios::binary);
  for(int i = 0; i <= lastStepInt; i++) {
  // Output this thread's data at each time step
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


  double gsize = GAMMA_ARRAY_SIZE, tsize = TAU_ARRAY_SIZE, csize = COAST_ARRAY_SIZE;
  output.open(outputPath + "finalOptimization-" + std::to_string(seed)+".bin", std::ios::binary);
  // output.open ("finalOptimization-"+std::to_string(static_cast<int>(seed))+"-"+std::to_string(threadRank)+".bin", std::ios::binary);

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

  // Optimized variables
  for (int j = 0; j < OPTIM_VARS; j++) {
    output.write((char*)&x[j], sizeof (double));
  }
  
  // Number of steps taken in final RK calculation
  output.write((char*)&lastStep, sizeof (double));

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
      terminalDisplay(adults[0], cConstants->missionObjectives);
  }

  //display to the terminal the best individual based on rankDistance
  std::cout << "\nBest Rank Distance Individual:";
  std::sort(adults.begin(), adults.end(), rankDistanceSort);
  terminalDisplay(adults[0], cConstants->missionObjectives);

  //Display number of errors
  std::cout << "\n# of errors this generation: " << numErrors << "\n";

  //Display number of duplicates
  std::cout << "\n# of duplicates this generation: " << numDuplicates << "\n";
  
  //display the oldest individual
  std::cout << "\nOldest age adult: " << generation - oldestBirthday << "\n";

  //Display the steps taken by the best individual
  std::cout << "\nBest step count: " << adults[0].stepCount << "\n";

  //Display the progress of the best rank distance individual 
  std::cout << "\nBest rank-distance adult progress: " << adults[0].progress << "\n\n";
}

// Utility function to display the currently best individual onto the terminal while the algorithm is still running
void terminalDisplay(const Adult& individual, const std::vector<objective> objectives) {
    
  //Print the parameters for each of the objectives for the passed in individual
  for (int i = 0; i < objectives.size(); i++) {
    //Print the name and value of the data of the objective
    std::cout << "\n\t" << objectives[i].name << ": " << individual.getParameters(objectives[i]);
  }

  std::cout << std::endl;
}