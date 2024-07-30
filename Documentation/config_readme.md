<h1> Config File Specifications/Information </h1>
<i>Last Updated: August 2022</i>

<h1>Document Links</h1>  

- [General Info](#general-info)
- [Values from genetic config](#values-from-genetic-config)
  - [Basic Variables](#basic-variables)
  - [GPU/Algorithim Variables](#gpualgorithim-variables)
  - [Display Variables](#display-variables)
  - [Input Parameter Range and Mutation Variables](#input-parameter-range-and-mutation-variables)
  - [Other Variables](#other-variables)
- [Mission config Values](#mission-config-values)
- [Target config Values](#target-config-values)
- [Impact Position and Velocity Values](#impact-position-and-velocity-values)
- [Recommended Config Values](#recommended-config-values)

<br>

# General Info

<h2>File and variable format for config file</h2>  

- The config file allows empty rows and comments ("//" at start of comment line) for formatting the presentation of the contents. The file also allows in-line comments with the only requirement being a space after the value assignment
- The parsing process currently does not attempt any verification when assigning values to variables (lack of assignment nor duplications).  When reading values, it takes the variable assignment and (for appropriate variable types) uses the standard string to int/double function. Variable assignment does not currently handle arithmetic operations
- When reading the file, assumptions are made that the config file contains valid simple values for all variables and that there are no spaces until after the variable value assignment

<h2>Config file structure</h2>  

- Whenever the code runs, three config files are used to gather info 
- First, genetic.config is read. This file holds info on the details of how the genetic algorithim should run, such as how strong mutations are. 
- Second, mission.config is opened. It holds the goals (Ex: relative speed)  of the genetic algorithim and the file name for the mission target file (Ex: psyche).
- Finally, the file containing info on the target body for the run is read. The information should include final position and velocities for the origin, target, and any graviational assist body. It should also include info on the setup of the spacecraft for the mission (like dry and fuel mass).

<h2>The cudaConstants struct</h2>

- In the code, the structure that uses the config file is called <b>cudaConstants</b> and is only accessed when being constructed (therefore changing the config file during a run would have no impact).
- Contains an overloaded << operator for cudaConstants that outputs the object's contents with labelling/formatting for better readibility to print onto the terminal screen in the main() function.
- For changing what config files are used, the file addresses can be changed in the cudaConstants constructor found in config.cpp. Requires re-compiling of the code.
- The default addresses are relative, based off of the .exe file created in the Optimization folder (Ex: "../Config_Constants/<name\>.config").
- If the config file address is invalid, the program will output to terminal that this is the case.  
<br>

# Values from genetic config
## Basic Variables
| Variable Name              	| Data Type  	| Units 	| Usage                                                                                                                                                      	                    |   	|
|----------------------------	|------------	|-------	|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---	|
| time_seed                  	| int/string 	| None  	| Sets the seed used in random generation initialization and labeling file outputs, either specify a seed to use or place "NONE" for the seed to be time(0)                         |   	|
| max_generations               | int           | None      | Determines how many generations the genetic algorithim will run before "giving up" and trying again on a new run   |       |
| run_count                     | int           | None      | The number of runs that will be completed before ending the program   |       |
| carryover_individuals         | int           | None      | The number of individuals from the previous run to use as a basis for the starting parameters of the initial children of the current run   |       |
| algorithm_type                     | string           | None      | What genetic diversity algorithm will be used. Either "rank-rarity" or "rank-distance"   |       |

## GPU/Algorithim Variables
| Variable Name              	| Data Type  	| Units 	| Usage                                                                                                                                                      	                    |   	|
|----------------------------	|------------	|-------	|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---	|
| rk_tol                 	    | double     	| None  	| The relative/absolute (not sure which one it is) tolerance for the runge kutta algorithm |   	|
| doublePrecThresh              | double     	| None  	| The smallest error difference in the runge kutta algorithm allowed. Setting the value too small would result in differences between GPU and CPU runge-kutta due to the data types limits of precision |   	|
| max_numsteps                 	| int        	| None  	| The maximum number of steps in the runge kutta used in the GPU and in EarthInfo's reverse runge kutta method. Note that changing numsteps but keeping the same time_seed value can lead to different results due to slight variance in the results that changes the algorithm's path to finding a solution    |   	|
| min_numsteps                 	| int        	| None  	| The minimum number of steps in the runge kutta used in the GPU and in EarthInfo's reverse runge kutta method. Note that changing numsteps but keeping the same time_seed value can lead to different results due to slight variance in the results that changes the algorithm's path to finding a solution    |   	|
| coast_threshold             	| double     	| None  	| In a range from 0 to 1. 1 sets the spacecraft to coast at all times while 0 sets the spacecraft to always have thruster on    |   	|
| timeRes                    	| int        	| seconds   | The "gap" between each calculation for Earth's backward runge-kutta, for example 3600 sets every calculation to be 1 hour apart   |   	|
| maxSimNum                    	| int        	| None   | Determines the maximum number of times to simulate each child before stopping the simulation and returning an error   |   	|
| num_individuals           	| int        	| None  	| Sets the size of the population pool and the number of threads used (each individual is given a separate thread). Recommended not to change |   	|
| survivor_count *depreciated*               	| int        	| None  	| Number of individuals selected as "survivors" to produce new individuals in the next generation of the genetic algorithm. Every pair produces 6 new individuals, value must be even. Right now, it's automatically set to 1/4 of num_individuals   |   	|
| thread_block_size           	| int        	| None  	| Number of threads per block on the GPU being used. Recommended to not change  |   	|

## Display Variables
| Variable Name              	| Data Type  	| Units 	| Usage                                                                                                                                                      	                    |   	|
|----------------------------	|------------	|-------	|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---	|
| record_mode                 	| boolean       | None  	| If "true", sets program to output various files that describe the performance. Meant to be used in helping to verify/debug behavior |   	|
| write_freq                 	| int        	| None  	| Sets number of generations to process before writing the current generation's information onto the genPerformance file. 1 is to write every generation   |   	|
| all_write_freq                | int           | None      | Sets the number of generations to process before writing information on all individuals in a generation onto files, (Ex:200 creates a new file every 200 generations)  |   |
| disp_freq                  	| int        	| None  	| Sets number of generations to process before outputting generation information to the console terminal. 1 is to display output every generation   |   	|

## Input Parameter Range and Mutation Variables
| Variable Name              	| Data Type  	| Units 	| Usage                                                                                                                                                      	                    |   	|
|----------------------------	|------------	|-------	|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---	|
| gamma_random_start_range      | double     	| None      | The magnitude of the value range (+/- the start range) for each child's random gamma coefficient initial value    |   	|
| tau_random_start_range        | double     	| None      | The magnitude of the +/- value range for tau coefficient random initial values    |   	|
| coast_random_start_range      | double     	| None      | The magnitude of the +/- value range for coast coefficient random initial values  |   	|
| alpha_random_start_range      | double     	| Radians   | The magnitude of the +/- value range for alpha random initial values  |   	|
| beta_random_start_range       | double     	| Radians   | The magnitude of the positive only value range for beta random initial values |   	|
| zeta_random_start_range       | double     	| Radians   | The magnitude of the +/- value range for zeta random initial values   |   	|
| mutation_amplitude            | double        | None      | A blanket scaling factor for the max range of mutations. This scalar is applied to any mutated value. A value of 1 means that this value does not modify mutations  |   |
| default_mutation_chance       | double    	| None  	| The probability of a mutations occurring when generating a new individual. Checks the mutation_chance before setting a random gene to be mutated and continues checking to mutate other genes until the check fails |   	|
| gamma_mutate_scale           	| double     	| None  	| Affects the maximum mutation range for gamma values (maximum mutation for the corresponding parameter is [this scale] * annealing * mutationScale)    |   	|
| tau_mutate_scale           	| double     	| None  	| Affects the maximum mutation range for tau values (maximum mutation for the corresponding parameter is [this scale] * annealing * mutationScale) 	                                                |   	|
| coast_mutate_scale           	| double     	| None  	| Affects the maximum mutation range for coast values (maximum mutation for the corresponding parameter is [this scale] * annealing * mutationScale) 	                                            |   	|
| triptime_mutate_scale 	    | double     	| Years  	| Affects the maximum mutation range for triptime values (maximum mutation for the corresponding parameter is [this scale] * annealing * mutationScale). When commented out, this value is set to the difference between triptime_min and triptime_max 	                    |   	|
| zeta_mutate_scale          	| double     	| Radians  	| Affects the maximum mutation range for zeta values (maximum mutation for the corresponding parameter is [this scale] * annealing * mutationScale) |   	|
| beta_mutate_scale           	| double     	| Radians  	| Affects the maximum mutation range for beta values (maximum mutation for the corresponding parameter is [this scale] * annealing * mutationScale) |   	|
| alpha_mutate_scale           	| double     	| Radians  	| Affects the maximum mutation range for alpha values (maximum mutation for the corresponding parameter is [this scale] * annealing * mutationScale)    |   	|

## Other Variables
| Variable Name              	| Data Type  	| Units 	| Usage                                                                                                                                                      	                    |   	|
|----------------------------	|------------	|-------	|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---	|
| anneal_initial             	| double     	| None  	| The initial anneal value used. Anneal impacts the maximum possible mutation value when generating a new individual (does not impact probability)  |   	|
| anneal_final             	  | double     	| None  	| The lowest that the anneal can reduce to. |  	 |
| best_count                 	| int        	| None  	| How many individuals must have obtained a solution before ending the algorithm, also outputs the top number of individuals up to best_count   |   	|
<br>

# Mission config Values
| Variable Name              	| Data Type  	| Units 	| Usage |   	|
|----------------------------	|------------	|-------	|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---	|
| Planetary_Objectives    | String        | N/A           | See [mission_config_readme.md](mission_config_readme.md) to find a detailed explanation on how to set this variable. This variable spot contains the requested objectives for the genetic algorithim. |   |
| Destination           | String       	| None  	| File file location for the requested mission or target. The linked file should have final position and velocity values for the origin and target body as well as information about the spacecraft |     |
| sun_r_min                 	| double       	| AU      	| Determines how close the spacecraft can get to the sun before being rejected.   |      |
| triptime_max                  | double     	| Years   | The maximum triptime, not only impacts the random starting guesses but also in deriving the time range calculation on PlanetInfo and also assumed to be larger than triptime_min |   	|
| triptime_min                  | double     	| Years   | The minimum triptime, not only impacts the random starting guesses but also in deriving the time range calculation on PlanetInfo and also assumes that it is less than triptime_max |   	|
<br>

# Target config Values
| Variable Name              	| Data Type  	| Units 	| Usage |   	|
|----------------------------	|------------	|-------	|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---	|
| Planetary Body Final Values   | double        | N/A       | These variables store the position and velocity of the mission's relevent planetary bodies at the time of mission completion. See the [Impact Position and Velocity Values](#impact-position-and-velocity-values) table for detail.  |   |
| gravAssistDist                | double        | AU        | The minimum distance a spacecraft can get from the gravitational assist body. |   |
| orbitalRadius                 | double        | AU        | The desired orbital radius for the spacecraft around the target. <b>Note: Include this variable in the target config file when only intending to do orbital missions!</b>    |   |
| orbitalSpeed                  | double        | AU        | The desired orbital speed for the spacecraft around the target. <b>Note: Include this variable in the target config file when only intending to do orbital missions!</b>    |   |
| thruster_type                	| int        	| None  	| Determine what thruster is used, 0 for none and 1 for NEXT ion thruster   |   	|
| dry_mass                     	| double        | kg      	| Set the mass of the spacecraft without fuel, also used in determining wet_mass    |   	|
| fuel_mass                     | double        | kg      	| Sets the initial mass of fuel in the spacecraft, used in determining wet_mass |   	|
| wet_mass                     	| double        | kg      	| The total mass of the spacecraft with fuel, value is derived after reading the config file when dry_mass and fuel_mass have had their values read |   	|
| c3scale                     	| double     	| None 	    | The scalar multiplier used to adjust c3energy by percentages for test runs, current implementation assumes that c3scale is assigned a value before c3energy (in config file c3scale is set before c3energy) 	|   	|
| c3energy                     	| double     	| m<sup>2</sup>/s<sup>2</sup>  	| The specific energy of the spacecraft when leaving the sphere of influence of the earth-moon center of mass, determines the magnitude of the escape velocity that is stored in v_escape 	|   	|
| v_escape                     	| double     	| AU/s  	| The magnitude of the initial velocity of the spacecraft when leaving the sphere of influence of the earth-moon center of mass, not in config file but rather derived from c3energy    |   	|
<br>


# Impact Position and Velocity Values
| Variable Name              	| Data Type  	| Units 	| Usage                                                                                                                                                      	                    |   	|
|----------------------------	|------------	|-------	|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---	|
| r_fin_target           	        | double     	| AU      	| The radius position of the target at impact date, relative to the Sun bodycenter  |   	|
| theta_fin_target        	        | double     	| Radians  	| The theta angle position of the target at impact date, relative to the Sun bodycenter |   	|
| z_fin_target           	        | double     	| AU      	| The z (off-plane offset) position of the target at impact date, relative to the Sun bodycenter    |   	|
| vr_fin_target           	        | double     	| AU/s  	| The velocity of the radius component of the target at impact date, relative to Sun bodycenter |   	|
| vtheta_fin_target      	        | double     	| AU/s  	| The tangental velocity of the target at impact date   |   	|
| vz_fin_target           	        | double     	| AU/s  	| The velocity of the z component of the target at impact date, relative to the Sun bodycenter  |   	|
| r_fin_earth           	    | double     	| AU     	| The radius position of the earth-moon center of mass at impact date, relative to the Sun bodycenter   |   	|
| theta_fin_earth          	    | double     	| Radians  	| The theta angle position of the earth-moon center of mass at impact date, relative to the Sun bodycenter  |   	|
| z_fin_earth           	    | double     	| AU      	| The z (off-plane offset) position of the earth-moon center of mass at impact date, relative to the Sun and used to plot it's path backwards in time for launch positions of the spacecraft    |   	|
| vr_fin_earth           	    | double     	| AU/s  	| The velocity of the radius component of the earth at impact date  |   	|
| vtheta_fin_earth         	    | double     	| AU/s  	| The tangental velocity of the earth-moon center of mass at impact date    |   	|
| vz_fin_earth           	    | double     	| AU/s  	| The velocity of the z component of the earth-moon center of mass at impact date, relative to Sun bodycenter and used to plot it's path backwards in time for launch positions of the spacecraft|   	|
| r_fin_mars             	    | double     	| AU     	| The radius position of the mars center of mass at impact date, relative to the Sun bodycenter |   	|
| theta_fin_mars          	    | double     	| Radians  	| The theta angle position of Mars' center of mass at impact date, relative to the Sun bodycenter   |   	|
| z_fin_mars               	    | double     	| AU      	| The z (off-plane offset) position of Mars' center of mass at impact date, relative to the Sun and used to plot it's path backwards in time for launch positions of the spacecraft |   	|
| vr_fin_mars              	    | double     	| AU/s  	| The velocity of the radius component of Mars at impact date   |   	|
| vtheta_fin_mars         	    | double     	| AU/s  	| The tangental velocity of Mars' center of mass at impact date |   	|
| vz_fin_mars           	    | double     	| AU/s  	| The velocity of the z component of Mars' center of mass at impact date, relative to Sun bodycenter and used to plot it's path backwards in time for launch positions of the spacecraft    |   	|
<br>

# Recommended Config Values
A table to elaborate on why some variables are assigned certain values

| Variable Name              	| Recommended Value  	| Reasoning 	|   |
|-----------------------------	|-------	| ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	|---    |
| [base bennu.config](\Base_Target_Configs\base_bennu.config)    | N/A   | The base config file used for missions where bennu is the target  |   |
| [base didymos.config](\Base_Target_Configs\base_didymos.config)    | N/A   | The base config file used for missions where didymos is the target  |   |
| [base marsOdyssey.config](\Base_Target_Configs\base_marsOdyssey.config)    | N/A   | The base config file used for missions where mars (oddesy mission) is the target  |   |
| [base marsOrbit.config](\Base_Target_Configs\base_marsOrbit.config)    | N/A   | The base config file used for missions where mars (orbital mission) is the target  |   |
| [base psyche.config](\Base_Target_Configs\base_psyche.config)    | N/A   | The base config file used for missions where psyche is the target  |   |
| mutation_rate                 | 0.5-0.75       | From test runs involving different mutation rates for the impact mission, 0.5 showed to be most effective at converging on a solution.  This was for singlue mutations, while mutation rates for more genes occurring after only showing no serious change. Rendezvous missions seem to need more mutations, so current rate 0.75.                         |       |
| triptime_min                  | 0.0-3.5         | These values are going to be mission dependent. For example, the bennu mission is best done as a redezvous mission with a triptime_min of 1.0. The didymos mission is similar. The psyche mission should have a min of 3.5, and the mars orbital should start at 0.5.   |       |
| triptime_max                  | 1.5-6.0       | As stated for the min, the max is also mission dependent. The bennu mission should end at 2.0. The didymos should end at 1.5. Psyche is a bit different as some of their missions plan to take 6 years while some end at 5 years. The mars orbital should be 1.5.  |       |
| timeRes                       | 3600      | Equivalent to 1 hour, this resolution for deriving earth location elements to store is considered sufficient.  For triptimes that fall between two indexes, the position/velocity is interpolated by a weighted average. This may need to be lower for mission that actaully use Mars' gravity as it is a very precise calculation.   |       |
| max_generations               | 5001-10001     | With current status of convergence rates for the algorithm in finding a valid solution being in the low thousands range, 10001 generations is considered plenty of time for the algorithm to find a solution and if it does reach this point then it won't find a solution as the annealing would become too small to lead to notable change that leads to a solution. Most missions should converge within 5000 generations. |       |
| num_individuals               | 2880      | Population size is based on the number of available threads in the Tesla GPU, optimizing the rate a pool can be calculated.   |       |
| survivor_count                | 720       | The genetic algorithim is able to dynamically create children as needed. However, through testing we have found 720 individuals to be a goood middle-ground value.    |       |
| thread_block_size             | 32        | Based on the Tesla GPU.   |       |
| best_count                    | 1         | It was found that in a given run, the best individuals when convergence occurs are very similar to each other, so setting more than 1 does not lead to more different solutions in a single run.  |       |
| doublePrecThresh              | 1e-12     | The double data type experiences loss in precision for differences in two elements within runge-kutta that as a result could lead to differences between the CPU and GPU computations of runge-kutta.  This value is a starting point in finding the smallest threshold for that difference. |     |