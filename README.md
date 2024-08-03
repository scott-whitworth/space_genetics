# Space Genetics
Based on previous work 2019 - 2021 Sankaran / Griffith summer research  
2019: Mateo, Lauren, Ben  
2020: Matthew, Noah, Andrew, Jeremy  
2021: Cassie, Trevor  
2022: Caleb, Connor, Meg  
2023: Caleb, Sierra  
2024: Josh

Pulled from: https://github.com/scott-whitworth/space_genetics

<h1>Multi-Objective Optimization Project</h1>
<i>Last updated: August, 2024</i>

<h2>Project Background & Current Objective</h2>

NASA has many missions to planets and asteroids with a variety of objectives: Rendezvous missions to collect samples of the target; orbital missions to study the target’s atmosphere; and missions where the spacecraft will perform a kinetic impact to change the orbital trajectory of its target. This program employs a multi-objective genetic algorithm which takes advantage of Nvidia's CUDA platform to optimize parameters to guide a spacecraft to its target. At this stage of development, the focus is to update the genetic algorithm to improve convergence on past missions and to allow the program to model more real missions, including those that rely on gravity assists.

Previously, the code could find the optimal trajectories for a spacecraft with a specified thruster to redirect an asteroid (as in the DART mission) or to perform a soft landing on the surface of an asteroid (as with the OSIRIS-REx mission). By contrast, for the summer of 2022, the goal was to improve the convergence rate of the current code and to adapt it to handle other mission types and optimize it for additional objectives (e.g. minimized fuel use or reduced trip time). A rank-distance genetic algorithm was implemented. The summer of 2023 focused on adding orbital missions and gravity assists, which required a new rank-rarity-based genetic algorithm. During the fall of 2023 to the spring of 2024, the objective portion of the algorithm was changed from maximization or minimization based objectives to objectives based on target values to allow for wider objective types. During the summer of 2024, the code was refined and verified to converge reliably for 3 three objective missions.


  Parameters being optimizing are the following
  | Variable Name               | Units    	  | Description                                                                                                                                                |   	|
  |----------------------------	|------------	|------------------------------------------------------------------------------------------------------------------------------------------------------------|---	|
  | Trip Time                  	| Seconds (s) | How long of a duration the spacecraft's trip takes once it leaves Earth's sphere of influence to impact with its target. The impact date is predetermined by the mission but the launch date is flexible so this variable impacts the initial position and velocity of the craft when it leaves Earth's sphere of influence|   	|
  | Alpha launch angle          | Radians  	  | The in-plane angle of the spacecraft as it leaves Earth's sphere of influence                                                                              |   	|
  | Beta launch angle           | Radians  	  | The in-plane angle from tangent of Earth's sphere of influence that the spacecraft's intial velocity takes                                                 |   	|
  | Zeta launch angle           | Radians  	  | The out-of-plane angle of the spacecraft's initial velocity when leaving Earth's sphere of influence                                                       |   	|
  | Gamma coefficients          | None  	    | Used in a fourier series to calculate the gamma (in-plane angle in radians) angle at a given time for thruster behavior                                    |   	|
  | Tau coefficients            | None  	    | Used in a fourier series to calculate the tau (out-of-plane angle in radians) angle at a given time for thruster behavior                                  |   	|
  | Coast coefficients          | None  	    | Used in a fourier series that determines if the thruster is activated or not (coasting) at a given time by comparing it to the coast_threshold |   	|

<h2>Files & Folders Navigation Guide</h2>

There are many folders and files in this project, so here's a short rundown of folder contents to help navigate everything;

  - Cuda: Where the most recent optimization code that attempts to find a best trajectory can be found. Uses the CUDA platform for GPU parallel processing in order to speed up the genetic algorithm.
    * [Config_Constants](./Documentation/config_readme.md): Where the cudaConstants structure is defined, the config files for different missions are stored, and where you specify what objectives to optimize for. The cudaConstants structure stores const values gathered from the necessary config files that we may want to change for different runs of the program. 
      * constants.h: Stores constant properties, such as AU unit value optimized variable offsets for the array that stores the values, and the masses and radii of planets and other celestial bodies. These are constants that should not be readily changed.
      * [Config files named after missions](./Documentation/config_readme.md#mission-config-values) (e.g. bennu.config): These files hold information about where Earth, the target, and Mars are when the mission ends as well as specific information about the spacecraft and other mission specific variables.
      * config.cpp/.h: Turns all the config files into values accessible through cConstants. If you are adding a new variable to a config, you need to declare it in config.h and set it up in config.cpp.
      * genetic.config: Holds max_numsteps, the timeseed, the number of individuals, the ranges for the different parameters, and things that relate to the genetic part of the algorithm and are less mission-specific.
      * [mission.config](./Documentation/mission_config_readme.md): Allows you to set the target and objectives you want the mission to optimize for.
    * [Planet_calculations](./Documentation/planet_calculations_info.md): Code for calculating the planet conditions and defines the global pointer variable launchCon (planetInfo.h). Dependent on Motion_Eqns/elements.h, Thrust_Files/thruster.h, and Config_Constants/config.h.
    * [Genetic_Algorithm](./Documentation/genetic_algorithm_info.md): Defines individuals used in the genetic algorithm, crossover/mutation methods to generate new generations in a pool, and the two ranking methods.
    * Motion_Eqns: Defines the elements structure that is used to describe the position and velocity of an object in space (such as Earth). Currently they take into account the gravitation of Mars and the Sun. Dependent on the Thrust_Files and Config_Constants/config.h.
    * Optimization: Contains the main file (optimization.cu) that has main() and the chief optimization function.
    * Output_Funcs: Contains the functions for outputting simulation information to files and the terminal display.
    * Runge_Kutta: Holds the runge_kutta functions with versions for both CPU and GPU usage. Also defines the rkParameters structure that holds the information the genetic algorithm attempts to optimize.  Dependent on Thrust_Files, Config_Constants, etc.
    * Thrust_Files: Contains the code that describes the thruster and its behavior using a Fourier series to determine angles of direction. Dependent on Config_Constants/config.h
  - PostProcessing: Contains MATLAB files that take in output files from the Cuda program and display further results.
  
<h2> Running CUDA Code: </h2>

1. To Compile:  
  On WU System & Personal Computers:
  1. Open VSCode and open the Space_Genetics project folder
  2. Open the command prompt in the VScode terminal (ctrl + tilde). Make sure you're under the "Terminal" tab.
  3. Navigate to the Optimization folder (input and enter `cd Cuda` then `cd Optimization`).
  4. Enter the following (new Summer 2021):
      `nvcc -ccbin "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.29.30133\bin\Hostx64\x64\cl.exe" -o test optimization.cu` 
    - `test` will be the name of the executable. This could be set to anything. 
    - *If your computer can't find the file, then there is a problem with the path. Copy and paste the file path into a file explorer, and delete everything after `\MSVC\` then click through until you get to `\cl.exe`. Copy this path and use it to replace the old path.*
    - This code should generate a .exe file in the Optimization folder with the output name that was entered

  For Tesla Machine:			  
   - Do the same as above, but for step 3 use this command:
     `nvcc -ccbin "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\bin\cl.exe" optimization.cu -o <nameOfOutput>`

2. To Run in VScode:
    1. Open command prompt in the VScode terminal (ctrl + tilde).
    2. Navigate to the Optimization folder (input and enter `cd Cuda` then `cd Optimization`).
    3. Type `.\` and then the .exe file name from above (don't forget to add .exe to the end) and hit enter.
    4. The program will begin and should perform the following processes in order;
       1. Output the GPU device name and intial values read from genetic.config (found in the Config_Constants folder).
       2. Calculate the Planet data with a visible loading bar. The range of this data is based on triptime_min and triptime_max in the target named config file
       3. On the terminal, displays a "." for every generation calculated and sorted.
       4. Every disp_freq generation, display the current generation number (how many have been calculated up to this point minus 1) along with the best progress individuals for each objective, best overall progress individual in the pool, and other generation information.
       5. Every write_freq generation, record that generation's performance, and every all_write_freq generation, record all individuals in that generation.
       6. Once the best individual(s) in the pool have passed the tolerance, the algorithm has "succeeded" and will output files that describe that individual's parameters that can be used in PostProcessing to analyze the run.
       7. Perform steps 2-6 again with a different time_seed value (offset by 100) if the number of runs performed is less than the config defined run value and the program has not succeeded yet.
       8. When the program has succeeded, it is finished and will close.

3. Changing properties:  
    In Config_Constants, the different config files hold a set of variables that can be changed before running the .exe file.  Refer to [config_readme.md](./Documentation/config_readme.md) for specifics on how each variable impacts the program's behavior and format of the file.  The code does not need to be recompiled if changes are made to a config file.
    The parameter lengths for gamma, tau, and coast can be manipulated in constants.h in lines 33-35 with the offsets handled based on those size values. The code must be recompiled for these changes to be put into effect. Additionally, changing the objective a mission is optimizing for is done in [mission.config](./Documentation/mission_config_readme.md).

4. Output Resulting Files (### refers to the time_seed value used in the run)
  - All these files from the run are put in a folder corresponding to their seed and the number of times the seed has been run: `<###>_<run times>`.
  - genPerformance-###.csv: Every write_freq generations, records the highest progress individual to the file along with the best objective values of that generation. Additionally the average objective values are recorded. Note: the best algorithm individual has its rkParameters output as well. Other metrics that could be helpful for troubleshooting (like the generation's anneal, statistics on the distances and ages of individuals in the generation, their progress, and how many duplicates and errors caused by getting too close to Mars) are also output.
  - simpleGenPerformance-###.csv: This file is a simplified version of genPerformance. Every write_freq generations, it records the progress of the best algorithm individual and the values of each objective of the best individuals by the genetic algorithm sort and by each separate objective. Note: printing of this file is currently disabled.
  - marsCheckValues-###-gen#.csv: Recorded prior to starting the optimization of the genetic algorithm. Stores the positions and velocity vectors calculated at a number of steps (currently each timeRes long).
  - earthCheckValues-###-gen#.csv: Recorded prior to starting the optimizing genetic algorithm. Contains the position and velocity for every day of Earth. Note: lower resolution than what is calculated (which is every hour).  Used to verify that it matches with the JPL database.
  - ###-AllAdults-gen#.csv: Every all_write_freq generations, records a lot of statistics on every individual in a generation. It records their positions, all their rkParameters, their ages, birthdays, ranks, distances, average parent progress, progress, and their parentChildProgressRatio.
  - referencePoints-###.csv: Prints the coordinates of all of the created reference points for the run.
  - runReports.csv: Once a run has finished -no matter if it converged- it is recorded here. Its seed, final generation, and whether or not it converged are recorded on one row. It then outputs basic statistics on the top 3 rankDistance individuals (their progress, age, and their performance on each objective).
  - finalOptimization-###.bin and orbitalMotion-###.bin: These files are created for the best progress individual when the runs are finished. This allows PostProcessing to show the found trajectory.

<h2> Printing orbital plots </h2>

1. In order to run the filePlot Matlab code, you must be on a Whitworth machine that has Matlab installed.
2. Copy the finalOptimization-###.bin and orbitalMotion-###.bin files outputted by the Cuda application when the run is finished into same folder as the matlab code being used (Under the Matlab_PostProcessing folder).
3. Open the Matlab client and navigate to the directory where the matlab code and binary files are placed (Under the Matlab_PostProcessing folder).
4. Run the matlab script in the command window:
   - To view an individual's path, enter "filePlot(#)" with # as the UTC time seed for the run.
   - To compare different paths, replace "filePlot(#)" with "filePlotCompare(#,#)" with the # being two different time seeds.
5. Graphs will be generated that show the path in three dimensions, the coast behavior, etc. and can be exported as pngs.

<h2> Printing reference point plots </h2>

1. NEEDS WORK

<h2>NASA JPL Data for Impact Date Position & Velocity</h2>
Here is how the impact, rendezvous, and orbital data was obtained to be used in the asteroid config files.

1.  Navigate to https://ssd.jpl.nasa.gov/horizons.cgi, this database is used to retrieve final position and velocity vector components for Earth, Mars, and the astroid's barycenter relative to the Sun.
2.  The ephemeris type should be a vector table. In table settings, select Type 2. Set the coordinate origin to the Sun's center of mass. The current impact date is 09-30-2022 19:54:55 UTC, so set an adequate time window with a resolution of 1 minute. Set the units to a.u. for position and a.u./day for velocity. Then, select Generate Ephemeris.
3.  To obtain final Earth elements, change the target body to Earth-Moon barycenter and generate a new ephemeris. To obtain the final elements for Mars, do not use the barycenter for Mars, just Mars as the target.
4.  Using [impactParams.m](Matlab_PostProcessing/impactParams.m), DerivingImpactValues.xlsx, or some other equivalent method, convert the values to cylindrical coordinates with velocity values changed from AU/day to AU/s.


