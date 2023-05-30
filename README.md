# Space Genetics
Based on previous work 2019 - 2021 Sankaran / Griffith summer research  
2019: Mateo, Lauren, Ben  
2020: Matthew, Noah, Andrew, Jeremy  
2021: Cassie, Trevor

2022: Caleb, Connor, Meg

Pulled from: https://github.com/scott-whitworth/Didymos-optimization on 6/7/2021  

<h1>Multi-Objective Optimization Project</h1>
<i>Last updated: August, 2022</i>

<h2>Project Background & Current Objective</h2>

NASA has many missions to planets and asteroids with a variety of objectives: Rendezvous missions to collect samples of the target; orbital missions to study the target’s atmosphere; and missions where the spacecraft will perform a kinetic impact to change the orbital trajectory of its target. This program employs a multi-objective genetic algorithm which takes advantage of Nvidia's CUDA platform to optimize parameters to guide a spacecraft fitted with a NEXT ion thruster to its target. At this stage of development, the focus is to update the genetic algorithm to optimize for more objectives (not only position difference and speed difference) and to allow the program to model more real missions, including orbital missions and those that rely on gravity assists.

Previously, the code could find the optimal trajectories for a spacecraft with the NEXT ion thruster to redirect an asteroid (as in the DART mission) or to perform a soft landing on the surface of an asteroid (as with the OSIRIS-REx mission). By contrast, for the summer of 2022, the goal was to improve the convergence rate of the current code and to adapt it to handle other mission types and optimize for additional objectives (e.g. minimized fuel use or reduced trip time). 


  Parameters being optimizing are the following
  | Variable Name               | Units    	  | Description                                                                                                                                                |   	|
  |----------------------------	|------------	|------------------------------------------------------------------------------------------------------------------------------------------------------------|---	|
  | Trip Time                  	| Seconds (s) | How long of a duration the trip for the spacecraft takes from once it leaves Earth's sphere of influence to impact with its target, the impact date is predetermined by the mission but launch date is flexible and so this variable also impacts the initial position and velocity of the craft when it leaves Earth's sphere of influence|   	|
  | Alpha launch angle          | Radians  	  | The in-plane angle of the spacecraft as it leaves Earth's sphere of influence                                                                              |   	|
  | Beta launch angle           | Radians  	  | The in-plane angle from tangent of Earth's sphere of influence that the spacecraft's intial velocity takes                                                 |   	|
  | Zeta launch angle           | Radians  	  | The out-of-plane angle of the spacecraft's initial velocity when leaving Earth's sphere of influence                                                       |   	|
  | Gamma coefficients          | None  	    | Used in a fourier series to calculate the gamma (in-plane angle in radians) angle at a given time for thruster behavior                                    |   	|
  | Tau coefficients            | None  	    | Used in a fourier series to calculate the tau (out-of-plane angle in radians) angle at a given time for thruster behavior                                  |   	|
  | Coast coefficients          | None  	    | Used in a fourier series that determines if the thruster is activated or not (coasting) at a given time for thruster behavior by comparing to coast_threshold |   	|

<h2>Files & Folders Navigation Guide</h2>

There are many folders and files in this project, so here's a short rundown of folder contents to help navigate everything;

  - Cuda: Where the most recent optimization code that attempts to find a best trajectory can be found, uses the CUDA platform to use  GPU and genetic algorithm 
    * [Config_Constants](../space_genetics/Cuda/Config_Constants/config_readme.md): Where cudaConstants structure is defined and the config files for different missions and that allow you to optmize for different objectives are. cudaConstants handles storing const values that we may want to be able to change for different runs of the program. 
      * constants.h: Stores constant properties, such as AU unit value optimized variable offsets for the array that stores the values, and the masses and radii of planets and other celestial bodies. These are constants that should not be easily changed.
      * [Config files name after missions](../space_genetics/Cuda/Config_Constants/config_readme.md#mission-config-values) (e.g. bennu.config): These files hold information about where Earth, the target, and Mars are when the mission ends, as well as specific information about the spacecraft and other mission specific variables.
      * config.cpp/.h: Turns all the config files into values accessible through cConstants. If you are adding a new variable to a config, you need to declare it in config.h and set it up in config.cpp.
      * genetic.config: Holds max_numsteps, the timeseed, the number of individuals, the ranges for the different parameters, and things that relate to the genetic part of the algorithm and are less mission-specific
      * [mission.config](Cuda/Config_Constants/mission_config_readme.md): Allows you to set the objectives you want the mission to optimize for.
    * [Planet_calculations](../Cuda/Planet_calculations/planet_calculations_info.md): Code for calculating the planet conditions and defines the global pointer variable launchCon (planetInfo.h). Dependent on Motion_Eqns/elements.h, Thrust_Files/thruster.h, and Config_Constants/config.h.
    * [Genetic_Algorithm](../space_genetics/Cuda/Genetic_Algorithm/genetic_algorithm_info.md): Defines individuals used in the genetic algorithm and crossover/mutation methods to generate new generations in a pool.
    * Motion_Eqns: Defines elements structure that is used to describe the position and velocity of an object in space (such as Earth). Currently they take into account the gravitation of Mars, as well as that of the Sun. Dependent on Thrust_Files and Config_Constants/config.h.
    * Optimization: Contains main file (optimization.cu) that has main() and main optimize function that is used.
    * Output_Funcs: Contains the functions for outputting files and terminal display, some methods such as recordAllIndividuals() are currently unused
    * Runge_Kutta: Holds the runge_kutta functions with versions for both CPU and GPU usage. Also defines rkParameters structure that holds the information that the genetic algorithm attempts to optimize.  Dependent on Thrust_Files, Config_Constants, etc.
    * Thrust_Files: Contains the code that describes the thruster and its behavior using a Fourier series to determine angles of direction. Dependent on Config_Constants/config.h
  - PostProcessing: Contains MATLAB files that take in output files from the Cuda program to display results.
  - DerivingImpactValues.xlsx : An old (likely temporary) method of converting JPL data values from cartesian to cylindrical units by using spreadhsheet fields and computations.  The converted values are then copied and put into genetic.config to be read and used.  This does not directly impact the code. Now there is a MATLAB function we generally use to do these calculations.
  
<h2> Running CUDA Code: </h2>

On WU System & Personal Computers:
1. To Compile:
   1. Open VSCode and open the Space_Genetics project folder
   2. Open command prompt in the VScode terminal, make sure this is `Command Line`.
   3. Navigate to the Optimization folder (input and enter "cd Cuda" then "cd Optimization").
   4. Enter the following (new Summer 2021):
      `nvcc -ccbin "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.29.30133\bin\Hostx64\x64\cl.exe" -o test optimization.cu -arch=compute_50 -code=sm_50` 
    - `test` will be the name of the executable. This could be set to anything.  
    - The last flags were added to force GPU to use its own architecture. This is dependent on which machine you are running. Those flags are for the machines in 214/304 Lab
    - *If your computer can't find the file, then there is a problem with the path. Copy and paste the file path into a file explorer, and delete everything after `\MSVC\` then click through until you get to `\cl.exe`. Copy this path and use it to replace the old path.*
    - This code should add an .exe file with the output name that was entered

  For Tesla Machine:			  
   - Do the same as above, but for step 3 use this command:
     `nvcc -ccbin "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\bin\cl.exe" optimization.cu -o <nameOfOutput>`

2. To Run in VScode:
    1. Open command prompt in the VScode terminal.
    2. Navigate to the Optimization folder (input and enter "cd Cuda" then "cd Optimization").
    3. Type the .exe file name from above (don't forget to add .exe) and enter
    4. The program will begin and should show the following in order on the terminal;
       1. Outputs the GPU device name and intial values read from genetic.config that is in Config_Constants folder.
       2. Calculate the Planet data with a visible loading bar.  The range of this data is based on triptime_min and triptime_max in config file
       3. Outputs the number of threads and blocks that will be used in the optimize function and starts the algorithm.
       4. On the terminal, displays a "." for every generation calculated and sorted.  Every disp_freq generation it displays the current generation number (how many have been calculated up to this point minus 1) and best speed, position and progress individuals in the pool.
       5. Along with terminal display, there are several file outputs made during the program's run.
       6. Once the best individual in the pool has passed the tolerance, the algorithm has "succeeded" and will output files that describe that individual's parameters that can be used in PostProcessing to observe.
       7. Perform steps 2-6 again with different time_seed value if the number of runs performed is less than the run value in the config.
       8. The program is finished and will close.

3. Changing properties:
      In Config_Constants, the different config files hold a set of variables that can be changed before running the .exe file.  Refer to [config_readme.md](Cuda/Config_Constants/config_readme.md) for specifics on how each variable impacts the program's behavior and format of the file.  The code does not need to be recompiled if changes are made to a config file.
      The parameter lengths for gamma, tau, and coast can be manipulated in constants.h in lines 23-25 with the offsets handled based on those size values. The code must be recompiled for these changes to be put into effect. Additionally, changing the objective a mission is optimizing for is done in mission.config.

4. Output Resulting Files (### refers to time_seed value used in the run)
  - All these files are put in a file corresponding to their seed and the number of times the seed has been run.
  - genPerformance-###.csv: Every write_freq generations, the objectives being optimized for of the best individuals by rankDistance sort and of the best individual for each objective are recorded to this file. Additionally the average posDiff and speedDiff are recorded. The best rankDistance individual has many of its rkParameters output as well. Other metrics that could be helpful for troubleshooting like the generation's anneal, statistics on the distances and ages of individuals in the generation, their progress, and how many duplicates and errors caused by getting too close to Mars are also output.
  - simpleGenPerformance-###.csv: This file is a simplified version of genPerformance. Every write_freq generations, it records the progress of the best rankDistance individual and the values of each objective of the best individuals by rankDistance sort and by each separate objective.
  - marsCheckValues-###-gen#.csv: Recorded prior to starting the optimizing genetic algorithm, the positions and velocity vectors calculated at a number of steps (currently each timeRes long).
  - earthCheckValues-###-gen#.csv: Recorded prior to starting the optimizing genetic algorithm, contains the position and velocity for every day of Earth (lower resolution than what is calculated which is every hour).  Used to verify that it matches with the JPL database.
  - ###-AllAdults-gen#.csv: Every all_write_freq generations, records a lot of statistics on every individual in a generation. It records their positions, all their rkParameters, their ages, birthdays, ranks, distances, average parent progress, progress, and their parentChildProgressRatio.
  - runReports.csv: Once a run has finished, whether or not it converged, it is recorded here. Its seed, final generation, and whether or not it converged are recorded on one row, then it outputs basic statistics on the top 3 rankDistance individuals (their progress, age, and their performance on each objective).
  - finalOptimization-###.bin and orbitalMotion-###.bin : For the best individual that has reached a solution, files are made of this format when the algorithm is finished. This allows PostProcessing to show the trajectory found.

<h2> Running Matlab Code </h2>

1. In order to run the Matlab code, you must be on a Whitworth machine that has Matlab installed.
2. Copy the binary files outputted by the Cuda application generated when it finished into same folder as the matlab code being used.
3. Open Matlab client with directory to where the matlab code files and binary files are placed.
4. Run the matlab script in the command window:
   - To view an individual's path, enter "filePlot(#)" with # as the UTC time seed for the run.
   - To compare different paths, replace "filePlot(#)" with "filePlotCompare(#,#)" with the # being two different time seeds.
5. Graphs will be generated that show the path in a three dimensional graph, coast behavior, etc. that could be exported into png format.

<h2>NASA JPL Data for Impact Date Position & Velocity</h2>
Here is how the impact, rendezvous, and orbital data was obtained to be used in the asteroid config files.

1.  Navigate to https://ssd.jpl.nasa.gov/horizons.cgi, this database is used to retrieve final position and velocity vector components for Earth, Didymos, Mars, and Bennu's barycenter relative to the Sun.
2.  The ephemeris type should be a vector table. In table settings, select Type 2. Set the coordinate origin to the Sun's center of mass. The current impact date is 09-30-2022 19:54:55 UTC, so set an adequate time window with a resolution of 1 minute. Set the units to a.u. for position and a.u./day for velocity. Then, select Generate Ephemeris.
3.  To obtain final Earth elements, change the target body to Earth-Moon barycenter and generate a new ephemeris. To obtain the final elements for Mars, do not use the barycenter for Mars, just Mars as the target.
4.  Using impactParams.m, DerivingImpactValues.xlsx, or some other equivalent method, convert the values to cylindrical coordinates with velocity values changed from AU/day to AU/s.

<h2>Flowchart Overview of CUDA Code</h2>
<i>Last updated on August 5th, 2021</i>
