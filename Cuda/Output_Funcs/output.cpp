#include <string>
#include <iomanip>
#include "math.h"

// Utility function to display the currently best individual onto the terminal while the algorithm is still running
// input: Individual to be displayed (assumed to be the best individual of the pool) and the value for the current generation iterated
// output: onto the console termina, generation is displayed and best individual's posDiff, speedDiff, and cost values
void terminalDisplay(Individual& individual, unsigned int currentGeneration) {
    std::cout << "\nGeneration: " << currentGeneration << std::endl;
    std::cout << "Best individual:" << std::endl;
    std::cout << "\tposDiff: " << individual.posDiff << std::endl;
    std::cout << "\tspeedDiff: " << individual.speedDiff << std::endl;
    std::cout << "\tcost: "    << individual.cost << std::endl;
}

// mutateFile[time_seed].csv is given a header row, now ready to be used by recordMutateFile()
void setMutateFile(const cudaConstants* cConstants) { 
  //Set up file
  std::ofstream mutateFile;
  int seed = cConstants->time_seed;
  mutateFile.open("mutateFile-" + std::to_string(seed) + ".csv", std::ios_base::app);

  //Set Headers
  mutateFile << "gen, anneal, genesToMutate,";
  for (int i = 0; i < GAMMA_ARRAY_SIZE; i++) {
    mutateFile << "gamma" << i << ",";
  }
  for (int i = 0; i < TAU_ARRAY_SIZE; i++) {
    mutateFile << "tau" << i << ",";
  }
  for (int i = 0; i < COAST_ARRAY_SIZE; i++) {
    mutateFile << "coast" << i << ",";
  }
  mutateFile << "alpha,beta,zeta,tripTime\n";

  mutateFile.close();
}

void errorCheck(double *time, elements<double> *yp,  double *gamma,  double *tau, int & lastStep, double *accel, double *fuelSpent, const double & wetMass, double *work, double *dE, double *Etot_avg, const cudaConstants* config) {
  // Initializing storage at each time step
  double *mass = new double[lastStep];
  double *Etot = new double[lastStep];
  
  //Iterate over all time steps
  for (int i = 0; i < lastStep; i++) {
    // total mass at each time step
    mass[i] = wetMass - fuelSpent[i];
    // total mechanical energy (K + U) at each time step
    Etot[i] = mass[i] * ((pow(yp[i].vr,2) + pow(yp[i].vtheta,2) + pow(yp[i].vz,2))/ 2 - constG * massSun / yp[i].r) / pow(AU,2);
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

// file mutateFile-[time_seed].csv is appended a new row containing mutate information
void recordMutateFile(const cudaConstants * cConstants, double generation, double annealing, int numGenes, double recordLog[OPTIM_VARS]) {
  std::ofstream mutateFile;
  int seed = cConstants->time_seed;
  mutateFile.open("mutateFile-" + std::to_string(seed) + ".csv", std::ios_base::app);

  // Record generation, annealing, and number of genes that should be impacted
  mutateFile << generation << "," << annealing << "," << numGenes << ",";

  // Record gamma mutation values
  for (int i = GAMMA_OFFSET; i < (GAMMA_OFFSET + GAMMA_ARRAY_SIZE); i++) {
      mutateFile << recordLog[i] << ",";
  }
  // Record tau mutation values
  for (int i = TAU_OFFSET; i < (TAU_OFFSET + TAU_ARRAY_SIZE); i++) {
      mutateFile << recordLog[i] << ",";
  }
  // Record coast mutation values
  for (int i = COAST_OFFSET; i < (COAST_OFFSET + COAST_ARRAY_SIZE); i++) {
      mutateFile << recordLog[i] << ",";
  }
  // Record alpha, beta, zeta, tripTime
  mutateFile << recordLog[ALPHA_OFFSET] << "," << recordLog[BETA_OFFSET] << "," << recordLog[ZETA_OFFSET] << "," << recordLog[TRIPTIME_OFFSET] << ",";
  mutateFile << "\n";
  
  mutateFile.close();
}

// Main Output, final results of genetic algorithm
// input: x[] - array of OPTIM_VARS for a single individual
//        generation - generation num of individual       
//        cConstants - Access constants info such as target element, earth element, derive spaceCraft element, also other values such as rk_tol
//        best - To access the best individual (pool[0])
// output: file orbitalMotion-[time_seed].bin is created that holds spacecraft RK steps and error
//         file finalOptimization-[time_seed].bin is created that holds earth/ast/ and trajectory parameter values
void trajectoryPrint( double x[], int generation, const cudaConstants* cConstants, Individual best) {
  /*set the asteroid and inital conditions for the earth and spacecraft:
  constructor takes in radial position(au), angluar position(rad), axial position(au),
  radial velocity(au/s), tangential velocity(au/s), axial velocity(au/s)*/

  // setting landing conditions of the asteroid (Sept 30, 2022)
  elements<double> asteroid = elements<double>(cConstants->r_fin_ast, cConstants->theta_fin_ast, cConstants->z_fin_ast, cConstants->vr_fin_ast, cConstants->vtheta_fin_ast, cConstants->vz_fin_ast);

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
  double timeFinal=cConstants->orbitalPeriod; // Orbital period of asteroid(s)
  double deltaT; // time step
  int numSteps = 5000; // initial guess for the number of time steps, guess for the memory allocated 
  deltaT = (timeFinal-timeInitial) / cConstants->GuessMaxPossibleSteps; // initial guess for time step, small is preferable

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
  rk4sys(timeInitial, x[TRIPTIME_OFFSET] , times, spaceCraft, deltaT, yp, absTol, coeff, gamma, tau, lastStepInt, accel_output, fuelSpent, wetMass, cConstants);

  // store the number of steps as a double for binary output
  double lastStep = lastStepInt;

  // gets the final y values of the spacecrafts for the cost function.
  elements<double> yOut = yp[lastStepInt];

  // calculate the error in conservation of mechanical energy due to the thruster
  errorCheck(times, yp, gamma, tau, lastStepInt, accel_output, fuelSpent, wetMass, work, dE, Etot_avg, cConstants);


  


  // This function is used to compare the final best thread with other runs
  // append this thread's info to a csv file
  if (cConstants->record_mode == true) {
    // Record initial and final fuel masses along with tripTime and relative velocity at impact
    recordFuelOutput(cConstants, x, fuelSpent[lastStepInt], best);
    progressiveAnalysis(generation, lastStepInt, x, yOut, cConstants);
  }

  // binary outputs
  std::ofstream output;
  int seed = cConstants->time_seed;
  output.open("orbitalMotion-"+std::to_string(seed)+".bin", std::ios::binary);
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
  output.open("finalOptimization-"+std::to_string(seed)+".bin", std::ios::binary);
  // output.open ("finalOptimization-"+std::to_string(static_cast<int>(seed))+"-"+std::to_string(threadRank)+".bin", std::ios::binary);

  // Impact conditions
  output.write((char*)&cConstants->r_fin_ast, sizeof(double));
  output.write((char*)&cConstants->theta_fin_ast, sizeof(double));
  output.write((char*)&cConstants->z_fin_ast, sizeof(double));
  output.write((char*)&cConstants->vr_fin_ast, sizeof(double));
  output.write((char*)&cConstants->vtheta_fin_ast, sizeof(double));
  output.write((char*)&cConstants->vz_fin_ast, sizeof(double));

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

// Record progress of individual
// input: output - the output file stream being used
//        rank - the positional performance of the individual
//        ind - the individual object being recorded
//        config - cudaConstants object for accessing thruster_type information
// output: output file is appended information on rank, individual values/parameter information
void progressiveAnalysis(int generation, int numStep, double *start, elements<double> & yp, const cudaConstants *config) {
  int seed = config->time_seed;
  //Set up file
  std::ofstream output;
  output.open("progressiveAnalysis.csv", std::ios_base::app);
  output << "\ntime_seed,numStep,posDiff,speedDiff,tripTime,alpha,beta,zeta,";

  //Headers
  for (int i = 0; i < GAMMA_ARRAY_SIZE; i++) {
    output << "gamma"; 
    if (i == 0) {
      output << "_a" << i << ",";
    }
    else if (i % 2 == 0) {
      output << "_b" << i/2 << ",";
    }
    else {
      output << "_a" << (i+1)/2 << ",";
    }
  }
  for (int i = 0; i < TAU_ARRAY_SIZE; i++) {
    output << "tau"; 
    if (i == 0) {
      output << "_a" << i << ",";
    }
    else if (i % 2 == 0) {
      output << "_b" << i/2 << ",";
    }
    else {
      output << "_a" << (i+1)/2 << ",";
    }
  }
  for (int i = 0; i < COAST_ARRAY_SIZE; i++) {
    output << "psi"; 
    if (i == 0) {
      output << "_a" << i << ",";
    }
    else if (i % 2 == 0) {
      output << "_b" << i/2 << ",";
    }
    else {
      output << "_a" << (i+1)/2 << ",";
    }
  }
  output << "\n";  

  output << seed << ',' << numStep << ','; 
  output << sqrt(pow(config->r_fin_ast - yp.r, 2) + pow(config->r_fin_ast * config->theta_fin_ast - yp.r * fmod(yp.theta, 2 * M_PI), 2) + pow(config->z_fin_ast - yp.z, 2)) << ',';
  output << sqrt(pow(config->vr_fin_ast - yp.vr, 2) + pow(config->vtheta_fin_ast - yp.vtheta, 2) + pow(config->vz_fin_ast - yp.vz, 2)) << ',';
  output << start[TRIPTIME_OFFSET] << ',' << start[ALPHA_OFFSET] << ',' << start[BETA_OFFSET] << ',' << start[ZETA_OFFSET] << ',';

  for (int i = GAMMA_OFFSET; i < GAMMA_ARRAY_SIZE + GAMMA_OFFSET; i++) {
    output << start[i] << ","; 
  }
  for (int i = TAU_OFFSET; i < TAU_ARRAY_SIZE + TAU_OFFSET; i++) {
    output << start[i] << ","; 
  }
  for (int i = COAST_OFFSET; i < COAST_ARRAY_SIZE + COAST_OFFSET; i++) {
    output << start[i] << ","; 
  }

  output << std::endl;
  output.close();
}

// Initialize the .csv file with header row
// input: cConstants - to access time_seed for deriving file name conventions and also thruster type
// output: file genPerformanceT-[time_seed].csv, is appended with initial header row info
void initializeRecord(const cudaConstants * cConstants) {
  std::ofstream excelFile;
  int seed = cConstants->time_seed;
  // setting the numeric id tag as the randomization seed (when doing runs of various properties, suggested to add other values to differentiate)
  std::string fileId = std::to_string(seed);
  excelFile.open("genPerformance-" + fileId + ".csv", std::ios_base::app);

  excelFile << "gen,lowerPosDiff,bestSpeedDiff,bestCost,alpha,beta,zeta,tripTime,";
  
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

  excelFile << ",\n";
  excelFile.close();
}

// Take in the current state of the generation and appends to excel file, assumes initializeRecord() had already been called before (no need to output a header row)
void recordGenerationPerformance(const cudaConstants * cConstants, Individual * pool, double generation, double new_anneal, int poolSize) {
  std::ofstream excelFile;
  int seed = cConstants->time_seed;
  std::string fileId = std::to_string(seed);
  excelFile.open("genPerformance-" + fileId + ".csv", std::ios_base::app);
  excelFile << std::setprecision(20);
  excelFile << std::fixed;
  // Record best individuals best posDiff and speedDiff of this generation
  excelFile << generation << "," << pool[0].posDiff << ",";
  excelFile << pool[0].speedDiff << ",";
  excelFile << pool[0].cost << ",";

  // Record best individuals parameters
  excelFile << pool[0].startParams.alpha << ",";
  excelFile << pool[0].startParams.beta << ",";
  excelFile << pool[0].startParams.zeta << ",";
  excelFile << pool[0].startParams.tripTime << ",";

  for (int i = GAMMA_OFFSET; i < GAMMA_ARRAY_SIZE + GAMMA_OFFSET; i++) {
    excelFile << pool[0].startParams.coeff.gamma[i-GAMMA_OFFSET] << ","; 
  }
  for (int i = TAU_OFFSET; i < TAU_ARRAY_SIZE + TAU_OFFSET; i++) {
    excelFile << pool[0].startParams.coeff.tau[i-TAU_OFFSET] << ","; 
  }
  for (int i = COAST_OFFSET; i < COAST_ARRAY_SIZE + COAST_OFFSET; i++) {
    excelFile << pool[0].startParams.coeff.coast[i-COAST_OFFSET] << ","; 
  }

  //New anneal every gen
  excelFile << new_anneal << ",";
  excelFile << "\n"; // End of row
  excelFile.close();
  
}

// Takes in a pool and records the parameter info on all individuals, currently unused
// input: cConstants - to access time_seed in deriving file name
//        pool - holds all the individuals to be stored
//        poolSize - to use in iterating through the pool
//        generation - used in deriving file name
// output: file generation#[generation]-[time_seed].csv is created with each row holding parameter values of individuals
void recordAllIndividuals(const cudaConstants * cConstants, Individual * pool, int poolSize, int generation) {
  std::ofstream entirePool;
  int seed = cConstants->time_seed;
  entirePool.open("generation#" + std::to_string(generation) + "-" + std::to_string(seed) + ".csv");
  // Setup the header row
  entirePool << "position,alpha,beta,zeta,tripTime,";
  for (int i = 0; i < GAMMA_ARRAY_SIZE; i++) {
    entirePool << "gamma" << i << ",";
  }
  for (int i = 0; i < TAU_ARRAY_SIZE; i++) {
    entirePool << "tau" << i << ",";
  }
  for (int i = 0; i < COAST_ARRAY_SIZE; i++) {
    entirePool << "coast" << i << ",";
  }
  entirePool << "cost,";
  entirePool << '\n';

  entirePool << std::setprecision(20);
  entirePool << std::fixed;

  // Record all individuals in the pool
  for (int i = 0; i < poolSize; i++) {
    entirePool << i << ",";
    entirePool << pool[i].startParams.alpha << ",";
    entirePool << pool[i].startParams.beta << ",";
    entirePool << pool[i].startParams.zeta << ",";
    entirePool << pool[i].startParams.tripTime << ",";

    for (int j = 0; j < GAMMA_ARRAY_SIZE; j++) {
      entirePool << pool[i].startParams.coeff.gamma[j] << ",";
    }
    for (int j = 0; j < TAU_ARRAY_SIZE; j++) {
      entirePool << pool[i].startParams.coeff.tau[j] << ",";
    }
    for (int j = 0; j < COAST_ARRAY_SIZE; j++) {
      entirePool << pool[i].startParams.coeff.coast[j] << ",";
    }
    entirePool << pool[i].cost << ",";
    entirePool << "\n";
  }
  entirePool.close();
}

// Method for doing recording information at the end of the optimization process
// input: cConstants - access record_mode, if record_mode == true then call progressiveRecord method, also passed into writeTrajectoryToFile method as well as progressiveRecord
//        pool - To access the best individual (pool[0])
//        generation - to record the generation value 
//        thrust - passed into progressiveRecord and writeTrajectoryToFile
// output: writeTrajectoryToFile is called, if in record_mode then progressiveRecord is called as well
void finalRecord(const cudaConstants* cConstants, Individual * pool, int generation) {
  // To store parameter values and pass onto writeTrajectoryToFile
  double *start = new double[OPTIM_VARS];

  // Output the final best individual
  for (int j = 0; j < pool[0].startParams.coeff.gammaSize; j++) {
      start[GAMMA_OFFSET + j] = pool[0].startParams.coeff.gamma[j];
  }
  for (int j = 0; j < pool[0].startParams.coeff.tauSize; j++) {
      start[TAU_OFFSET + j] = pool[0].startParams.coeff.tau[j];
  }
  for (int j = 0; j < pool[0].startParams.coeff.coastSize; j++) {
      start[COAST_OFFSET + j] = pool[0].startParams.coeff.coast[j];
  }

  start[TRIPTIME_OFFSET] = pool[0].startParams.tripTime;
  start[ALPHA_OFFSET] = pool[0].startParams.alpha;
  start[BETA_OFFSET] = pool[0].startParams.beta;
  start[ZETA_OFFSET] = pool[0].startParams.zeta;

  // Test outputs
  std::cout << "Comparison\n";
  std::cout << "CUDA posDiff: " << pool[0].posDiff << std::endl;
  std::cout << "CUDA speedDiff: " << pool[0].speedDiff << std::endl;

  // Evaluate and print this solution's information to binary files
  trajectoryPrint(start, generation, cConstants, pool[0]);

  // cleaning up dynamic memory
  delete [] start;
}

// Stores information of launchCon of timeRes*24 resolution in EarthCheckValues-[time_seed].csv
void recordEarthData(const cudaConstants * cConstants) {
  double timeStamp = cConstants->triptime_min; // Start the timeStamp at the triptime min (closest to impact)
  // seed to hold time_seed value to identify the file
  int seed = cConstants->time_seed;
  // Open file
  std::ofstream earthValues;
  earthValues.open("EarthCheckValues-"+ std::to_string(seed) +".csv");
  // Set header row for the table to record values, with timeStamp
  earthValues << "TimeStamp, Radius, Theta, Z, vRadius, vTheta, vZ\n";
  // Fill in rows for elements of launchCon across the time range going backwards until triptime max
  while (timeStamp < cConstants->triptime_max) {
      earthValues << timeStamp << "," << launchCon->getCondition(timeStamp);
      timeStamp += cConstants->timeRes*24; // Increment to next day as timeRes is assumed to be set to every hour in this output
  }
  
  // Done recording earth calculations, close file
  earthValues.close();
}

// Record initial and final fuel masses along with tripTime and relative velocity at impact
// input: cConstants - access config constants
//        solution - best individual parameters from the final pool
//        fuelSpent - total fuel spent
//        best - To access the best individual (pool[0])
// output: fuelOutput.csv - output file holding fuel consumption and impact data
void recordFuelOutput(const cudaConstants* cConstants, double solution[], double fuelSpent, Individual best) {

  std::ofstream excelFile;
  excelFile.open("fuelOutput.csv", std::ios_base::app);
  
  excelFile << "\nSeed,Initial fuel (kg),Spent fuel (kg),Trip time (days),Impact speed (m/s)\n";

  excelFile << cConstants->time_seed << ",";
  excelFile << cConstants->fuel_mass << ",";
  excelFile << fuelSpent << ",";
  excelFile << solution[TRIPTIME_OFFSET]/(3600*24) << ",";
  excelFile << best.speedDiff*AU << "\n";

  excelFile.close();
}
