<h1>Code Overview (CS)</h1>
<i>This document is intended to explain the basic purposes behind the functions and classes used throughout the code, specifically though a CS perspective</i>
<br>
<i>Last Updated: Aug 2024</i>
<br>
</br>

# Table of Contents
- [High Level Overview](#high-level-overview)
- [Files](#files)
  - [Optimization.cu](#optimizationcu)
    - [main](#main)
    - [optimize](#optimize)
    - [checkTolerance](#checktolerance)
    - [calculateGenerationValues](#calculategenerationvalues)
  - [config.cpp/h](#configcpph)
  - [constants.h](#constantsh)
  - [Anneal.cpp/h](#annealcpph)
    - [changeAnneal](#changeanneal)
  - [ga_crossover.cpp/h](#ga_crossovercpph)
    - [maskValue enum](#maskvalue-enum)
    - [crossOver_randHalf](#crossover_randhalf)
    - [crossOver_oneParent](#crossover_oneparent)
    - [crossOver_wholeRandom](#crossover_wholerandom)
    - [crossOver_bundleVars](#crossover_bundlevars)
    - [crossOver_average](#crossover_average)
    - [crossOver_averageRatio](#crossover_averageratio)
    - [flipMask](#flipmask)
    - [copyMask](#copymask)
    - [getRand](#getrand)
    - [mutateMask](#mutatemask)
    - [mutate](#mutate)
    - [generateNewChild](#generatenewchild)
    - [generateChildrenPair](#generatechildrenpair)
    - [generateChildrenFromCrossover](#generatechildrenfromcrossover)
  - [geneticAlgorithm.cpp/h](#geneticalgorithmcpph)
    - [newGeneration](#newgeneration)
    - [fillParents](#fillparents)
    - [makeChildren](#makechildren)
    - [convertToAdults](#converttoadults)
    - [createFirstGeneration](#createfirstgeneration)
    - [firstGeneration](#firstgeneration)
    - [preparePotentialParents](#preparepotentialparents)
    - [addToAllAdults](#addtoalladults)
  - [Sort.cpp/h](#sortcpph)
    - [giveRank](#giverank)
    - [giveDistance](#givedistance)
    - [parameterSort](#parametersort)
  - [planetInfo.cpp/h](#planetinfocpph)
    - [PlanetInfo](#planetinfo)
    - [getCondition/getConditionDev](#getconditiongetconditiondev)
    - [calcIndex](#calcindex)
    - [calc_time](#calc_time)
    - [getTolData](#gettoldata)
    - [getAllPositions](#getallpositions)
    - [getPlanetPosition](#getplanetposition)
    - [interpolate](#interpolate)
    - [~PlanetInfo](#planetinfo-1)
    - [planetInitial_incremental](#planetinitial_incremental)
    - [getPlanetSize](#getplanetsize)
  - [ReferencePoints.cpp/h](#referencepointscpph)
    - [giveRarity](#giverarity)
    - [calcNormVector](#calcnormvector)
    - [calculateRelCost](#calculaterelcost)
    - [findPointDist](#findpointdist)
    - [findAssociatedPoints](#findassociatedpoints)
    - [calculateRarity](#calculaterarity)
    - [assignReservedRarity](#assignreservedrarity)
  - [motion_equations.cpp/h](#motion_equationscpph)
    - [calc_k and calc_kPlanet](#calc_k-and-calc_kplanet)
    - [calcRate_r](#calcrate_r)
    - [calcRate_theta](#calcrate_theta)
    - [calcRate_z](#calcrate_z)
    - [calcRate_vr](#calcrate_vr)
    - [calcRate_vtheta](#calcrate_vtheta)
    - [calcRate_vz](#calcrate_vz)
    - [calcRate_vrPlanet](#calcrate_vrplanet)
    - [calcRate_vthetaPlanet](#calcrate_vthetaplanet)
    - [calcRate_vzPlanet](#calcrate_vzplanet)
  - [output.cpp/h](#outputcpph)
    - [printBestAdults](#printbestadults)
    - [terminalDisplay](#terminaldisplay)
  - [runge_kutta.cpp/h](#runge_kuttacpph)
    - [rk4sys](#rk4sys)
    - [rk4Reverse](#rk4reverse)
    - [rkCalc](#rkcalc)
    - [rkCalcPlanet](#rkcalcplanet)
    - [calc_scalingFactor](#calc_scalingfactor)
    - [pmLimitCheck](#pmlimitcheck)
  - [runge_kuttaCUDA.cu/cuh](#runge_kuttacudacucuh)
    - [callRK](#callrk)
    - [rk4CUDASim](#rk4cudasim)
  - [calcFourier.cpp/h](#calcfouriercpph)
    - [calc_Series](#calc_series)
    - [calc_gamma](#calc_gamma)
    - [calc_tau](#calc_tau)
    - [calc_coast](#calc_coast)
- [Classes](#classes)
  - [Child](#child)
    - [Child_Members](#child_members)
    - [Constructors](#constructors)
    - [getParameters](#getparameters)
    - [getPosDiff](#getposdiff)
    - [getSpeedDiff](#getspeeddiff)
    - [getHorzVelDiff](#gethorzveldiff)
    - [getVertVelDiff](#getvertveldiff)
    - [getOrbitPosDiff](#getorbitposdiff)
    - [getOrbitSpeedDiff](#getorbitspeeddiff)
    - [getProgress](#getprogress)
  - [Adult](#adult)
    - [Adult_Members](#adult_members)
    - [Constructors](#constructors-1)
    - [bestProgress](#bestprogress)
    - [worstProgress](#worstprogress)
    - [lowAssocPoint](#lowassocpoint)
    - [dominationCheck](#dominationcheck)
    - [rankSort](#ranksort)
    - [rankDistanceSort](#rankdistancesort)
    - [rankRaritySort](#rankraritysort)
    - [duplicateCheck](#duplicatecheck)
    - [findDuplicates](#findduplicates)
  - [GPUMem](#gpumem)
    - [Members](#members)
    - [Initialize](#initialize)
    - [Free](#free)
  - [PlanetInfo(Class)](#planetinfoclass)
  - [Elements](#elements)
    - [Members](#members-1)
    - [Constructors](#constructors-2)
    - [operator+](#operator)
    - [operator-](#operator-1)
    - [operator*](#operator-2)
    - [operator/](#operator-3)
    - [compare](#compare)
    - [operator<<](#operator-4)
  - [rkParameters](#rkparameters)
    - [Description](#description-1)
    - [Members](#members-2)
    - [Constructors](#constructors-3)
    - [randomParameters](#randomparameters)
    - [compare](#compare-1)
    - [parametersRK4Simple](#parametersrk4simple)
  - [objective](#objective)
    - [Description](#description-1)
    - [Members](#members-3)
    - [Constructors](#constructors-4)
  - [ReferencePoints (class)](#referencepoints-class)
    - [Description](#description-2)
    - [Members](#members-4)
    - [Constructors](#constructors-5)
    - [addPoint](#addpoint)
  - [output (class)](#output-class)
    - [Members](#members-5)
    - [Constructor](#constructor)
    - [printGeneration](#printgeneration)
    - [printFinalGen](#printfinalgen)
    - [initializeGenPerformance](#initializegenperformance)
    - [initializeSimpleGenPerformance](#initializesimplegenperformance)
    - [prepareFolders](#preparefolders)
    - [recordGenerationPerformance](#recordgenerationperformance)
    - [recordGenSimple](#recordgensimple)
    - [recordAllIndividuals](#recordallindividuals)
    - [reportRun](#reportrun)
    - [trajectoryPrint](#trajectoryprint)
    - [errorCheck](#errorcheck)
    - [recordEarthData](#recordearthdata)
    - [recordMarsData](#recordmarsdata)
    - [recordReferencePoints](#recordreferencepoints)
  - [coefficients](#coefficients)
    - [Members](#members-6)
    - [initCoefficient](#initcoefficient)
    - [operator<<](#operator-4)



<br>
</br>

# High Level Overview
Here is a basic summary of how the genetic algorithm operates:  
1. The program will start in optimization.cu which first takes in the config values and then loops through each run by calling optimize()
2. In optimize(), the program creates a config defined number of individuals and simulates them (finding how well they performed) to fill a vector called `oldAdults`
3. The code then goes through a loop (the genetic algorithm) until an individual either converges or hits the maximum defined number of generations. In this loop:
    1. The algorithm calls newGeneration() which takes the individuals in `oldAdults` and fills up `newAdults` with oldAdult offspring that have gone through crossover and mutation
    2. The algorithm calls preparePotentialParents() which takes both `oldAdults` and `newAdults`, combines the error free individuals into a vector called `allAdults`, ranks these adults, and moves a config defined number of the highest ranking adults into the (now empty) `oldAdults` vector.
    3. One generation is now finished, so print a dot `.`
    4. Change the anneal based on the current best progress (which will affect the next batch of offspring's mutation range)
    5. Output statistics for this generation
4. Once the genetic algorithm has ended, print out the final generation data.

<br>
</br>

# Files

## Optimization.cu
### main
  - Function which starts the algorithm
  - Gathers CUDA device properties, reads the inputted config variables, calculates the possible positions for the Earth and Mars, calculates the reference points, and sets the run's base rng seed 
  - Runs a for loop which will call the [optimize()](#optimize) function for the number of runs specified in the config files
    - Recalculates the potential positions for planetary bodies each time it loops

### optimize
  - This function facilitates finding the optimal path for the input mission
  - It initially creates the 0th generation, made up of adults that are either imported from a previous run (See the carryover_individuals variable in [genetic.config](../Cuda/Config_Constants/genetic.config)) or are randomly generated
  - It then starts a do while loop which ends when either the generation limit has been reached or it has found a converged result
    - The first step in the loop is to generate new childeren from existing adults
    - Then, each individual is assigned a rank and a distance or rarity which is determined by being compared to all other individuals
    - Following this, generation-wide values are calculated for this generation
      - These values are primarily min/max/avg of outputs from the individuals' simulations
    - Then, the function will calculate the anneal for the next generation
    - Output files are then created or appended to with up-to-date information
    - Finally, there is a check to see if the best individual(s) have converged
  - Once best individual(s) number of individuals has converged or the generation limit has been reached, record final information

### checkTolerance
  - Checks the top x number of individuals to see if they have found an optimized answer
    - x (best_count) is defined in [genetic.config](../Cuda/Config_Constants/genetic.config)
  - It dynamically checks parameters depending on the number of input objectives
  - It returns true if the top x individuals have found an answer

###  calculateGenerationValues
  - Calculates generation-wide values which are passed to functions which write output files
  - List of values calculated:
    - Average, minimum, and maximum distance values
    - Average age, average birthday, and oldest birthday 
    - Average and best values for each objective
    - Average, min, and max steps
    - The number of reference points used

<br>

## config.cpp/h
- Sets us the constants that will be used throughout the code
- See [config_readme.md](config_readme.md) for more information on this section and for information on the config files themselves
<br>

## constants.h 
- Stores constants to be used throughout the code, along with enumerations to stores various statuses.
<br>


## Anneal.cpp/h
### changeAnneal
  - Calculates the anneal for the next generation
    - Anneal affects the range of the mutations for children parameters,. As the anneal gets smaller, mutation varies less which allows for smaller changes to input parameters when close to solving the objectives
  - Current implementation relies on a step method 
    - Each step sets the current anneal to a multiplier of a initial anneal value
    - The step is chosen based on the progress of the best individual
    - The higher the progress, the lower the anneal is set to
    - When the progress of the best individual gets high enough, the anneal calculation method switches to an exponential decay function
  - The base current anneal set by the steps is modified by a sinusoidal-based multiplier
    - As of writing, the current anneal is modified by up to +/- 70% of the base value of the current anneal
      - the modifier value is determined by a sinusoidal curve based on the current generation
  - Finally, it will make sure the new anneal never dips below the final anneal value

<br>

## ga_crossover.cpp/h
### maskValue enum
  - Used to generate masks for the crossover portion of generating children
  - Determines if a child's input parameter comes from parent 1 (1), parent 2 (2), the average of the parents (3), or a weighted average of the parents (4)
  - An array with the size of the number of input parameters is filled with instances of this enum; this is used to determine the origin of each of a child's input parameters

### crossOver_randHalf
  <i> Currently not in use</i>
  - Sets the the mask's indexes to be either parent 1 or parent 2
  - In this function, a cross index is randomly generated
  - Any index before the cross index is assigned to parent 1 and any index after is assigned to parent 2

### crossOver_oneParent
<i> Currently not in use</i>
  - Assigns all indexes of the mask array to parent 1

### crossOver_wholeRandom
  - Randomly assigns all indexes of the mask array to either parent 1 or parent 2

### crossOver_bundleVars
  - Randomly assigns the indexes of the mask array to either parent 1 or parent 2, but keeps the fourier series variables together.
  - This differs from crossOver_wholeRandom because the variables for all the gamma, tau, or coast parameters are just from one parent
    - Eg. All the gammas may be parent 1, the taus may be parent 2, and the coasts could all be parent 1

### crossOver_average
  - Assigns all indexes of the mask array to average

### crossOver_averageRatio
  - Assigns all indexes of the mask array to average ratio

### flipMask
  - Changes the enum value in each index of the mask array to its counterpart value
    - Swaps parent1 to parent2 and vice versa
    - Swaps average to average ratio and vice versa

### copyMask
<i> Currently not in use</i>
  - Takes in two crossover mask arrays and copies all of the enum values in the first/maskIn array to the second/maskOut array


### getRand
  - Generates a random number between the positive and negative value of an inputted maximum value
    - For example, if the argued maximum is 10, it will return a random value between 10 and -10
  - Used for generating input parameter mutations

### mutateMask
  - Fills in a mutation mask for a new child, which is an array of booleans that determines which parameters to mutate (when set to true)
  - The mask, with the size of the number of input parameters, is initially marked false for all input parameters
  - The function generates a random number between 0 and 1 and will mutate a parameter if the number is smaller than the mutation rate defined in genetic.config
    - If the number is larger, the function will stop marking parameters for mutation
  - If the number is smaller, the function decides a parameter needs to be mutated, and will randomly mark a new parameter to mutate by setting the corresponding index on the mutation mask to true.
    - It will keep randomly choosing parameters to mutate until the random number between 0-1 is greater than the mutation rate or all parameters are marked for mutation
  - Function used only by the [mutate()](#mutate) function

### mutate
  - Mutates children newly created from crossover
  - Uses a mutation mask that's altered in [mutateMask()](#mutatemask) to determine which of the child's parameters to mutate
  - To mutate a parameter it uses [getRand()](#getrand) and passes in the current [anneal](#annealcpph) which is multiplied by the input and global parameter mutate scales from [genetic.config](../Cuda/Config_Constants/genetic.config) to modify input parameters marked for mutation by a random +/- value
    - The anneal multiplier means the mutations will get smaller as the best individuals get closer to solving the objectives

### generateNewChild
  - Generates a set of input parameters for a new child based on a crossOver mask defined in [generateChildrenFromCrossover](#generatechildrenfromcrossover)
  - The parameters come from the two adults passed in
  - The mask is an array with indexes corresponding to each input parameter. Each index either has enum values of Partner1, Partner2, Avg, or Avg_Ratio, which will determine the origin of each of the child's input parameters
  - If the mask index for a parameter has Parent1 as the value, the parameter will be pulled from the first adult passed in
  - If the mask index for a parameter has Parent2 as the value, the parameter will be pulled from the second adult passed in
  - If the mask index for a parameter has Avg as the value, the parameter will be set as the 50/50 average of the parents' parameters
  - Finally, if the mask index for a parameter has Avg_ratio as the value, the parameter will be set as a weighted average of the parents' parameters
    - The weight of the average is randomly generated at the beginning of the function
  - After a set of parameters is generated, it will be mutated using the [mutate()](#mutate) function

### generateChildrenPair 
  - Creates up to two new children using generateNewChild and appends them to an array with new children
  - Flips the crossOver mask with [flipMask()](#flipmask) after generating the first child, so two unique children should be generated
  - Function will exit after generating the first child if the limit for children generation has been reached

### generateChildrenFromCrossover
  - Generates a number of children by the genetic crossover method, using a vector with [maskValue](#maskvalue-enum)s
    - The number of children generated comes from an inputted value 
  - Creates parent pairs by iterating through a shuffled list of the potential parents 
    - Note: The potential parents themselves are not shuffled, rather a vector with indexes corresponding to each future parent is shuffled
  - As of writing, it generates 6 children from three crossover methods
    - The used methods are crossOver_wholeRandom, crossOver_average, and crossOver_bundleVars
    - The same two parents generate all 6 new children
    - The number of children a pair of parents generates is flexible
  - If all parents have created children and children still need to be generated, the potential parents will be reshuffled, allowing for new parent pairs

<br>

## geneticAlgorithm.cpp/h
### newGeneration
  - Calls all functions necessary for generating a new set of children which will fill the newAdults vector
  - Functions called are:
    1. [fillParents()](#fillparents)- this selects the parents to be used when generating children
    2. [makeChildren()](#makechildren)- this creates the new set of children from the selected parents
    3. [callRK](#callrk)- this passes the children along to be simulated
    4. [convertToAdults](#converttoadults)- this fills the newAdults vector with the simulated children

### fillParents
  - Selects the top adults from the previous generations to be parents
    - The number of parents considered comes from the survivor_count variables in genetic.config
  - If it is generation 0, parents are selected from the first adults who were randomly generated
  - Potential parents will only include adults with error statues of valid or duplicate

### makeChildren
  - Currently, it only calls the [generatingChildrenFromCrossover](#generatechildrenfromcrossover) function, making sure the number of children generated equals num_individuals, which is defined in genetic.config

### convertToAdults
  - This first completes final output calculations and error checks for the new children
    - If it detects a child has nans for final position or velocity, it will assign an error to the child
    - It will then calculate outputs, some of which need their own dedicated [functions](#getposdiff)
  - Then, all new children are pushed to the newAdults vector (which is emptied at the beginning of the function)

### createFirstGeneration
  - This function creates a number of new children which are either randomly generated or pulled from a document and passes them to the [firstGeneration function](#firstgeneration)
    - If carryover_individuals is not set to 0 in [genetic.config](./config_readme.md/#values-from-genetic-config) and the program is able to open the allAdults file printed at the max generation from the file of the current time seed minus 100, the individual info will be taken from the opened allAdults file to make the initial children

### firstGeneration
  - Passes the children generated from createFirstGeneration to the simulation using [callRK](#callrk)
  - [Converts](#converttoadults) the simulated initial children into adults and pushes them into the oldAdults vector

### preparePotentialParents
  - Function facilitates the calculation of rank and distance (or rarity) for all individuals
    - This allows for the right parents to be chosen to create offspring for the next generation and for the best individual to be printed and checked for convergence
  - It will first determine which individuals are duplicates using [findDuplicates](#findduplicates)
  - Then it will add all valid adults to the allAdults vector using [addToAllAdults](#addtoalladults)
    - Note: the change which will get rid of or keep duplicates is within this function
  - Then it will calculate the rank and distance (or rarity) for all individuals in allAdults using the respective functions ([giveRank](#giverank) and [giveDistance](#givedistance) (for rank distance) ). If the genetic algorithm is set to rank rarity, instead of giveDistance(), it will call [giveRarity()](#giverarity) and [assignReservedRarity](#assignreservedrarity).
  - Finally, it will replace the individuals in the oldAdults vector with the best performing (rank distance or rank rarity) individuals from allAdults up to the config defined [num_individuals](./config_readme.md/#gpualgorithim-variables)

### addToAllAdults
  - Goes through all individuals in the inputted adult vector and checks their error statuses
    - If an adult's error status is not valid or duplicate, an error tracking variable is incremented. These are used for reporting
    - If an adult's error status is duplicate, one is added to a duplicate tracker used for reporting
      - Note: if you want to keep duplicates when they are detected, uncomment the line where duplicate adults are pushed to allAdults
    - If an adult's error status is valid, it is pushed to the allAdults vector

<br>

## sort.cpp/h
### giveRank
  - Assigns a rank to all valid (or duplicate) individuals
  - A rank can be thought of as tiers, where multiple individuals can occupy the same rank
    - However, when solving for only one parameter, rank works with each unique individual occupying a unique rank
  - Assigning a rank occurs in two stages:
  1. Each individual is compared against every other individual passed in and there is a check to see if one individual [dominates](#dominationcheck) the other
      - See the [section](#dominationcheck) on the adult class to learn about how an individual dominates another
      - Each individual is tracked for how many times they were dominated and what other individuals they have dominated
      - At the end of the comparisons, any individual who was never dominated is added to the "front" vector and given a rank of 1 (the best rank)
  1. Each individual is assigned a rank depending on how many times they were dominated
       - This process utilizes a while loop which runs until all individuals have a rank; the loop is:
         1. Go through all individuals in the front vector and pull up the individuals they dominated
         2. For each of the dominated individuals, subtract one from their counter for how many times they were dominated
         3. If an individual's domination counter reaches 0, add them to the "newFront" vector and assign them their rank (rank = times though this loop + 1)
         4. Once the front vector has been iterated through, replace the front vector with the newFront vector and clear the newFront vector
            - If the newFront is empty during this step, end the loop instead of repeating

### giveDistance
 - Assigns a distance to all valid (or duplicate) individuals. (nonvalid results get a distance of 0)
 - Distance measures how diverse an individual is compared to its neighbors
   - It is calculated by dividing the objective difference (between the parameter and target) of the individuals immediately better and worse than the focus individual by a normalization value
     - The normalization value is always the largest value for the objective difference within the population of individuals
 - First counts the number of valid or duplicate individuals so only they will be compared
 - Makes sure all distances are reset to 0 before calculating new distances
 - Note: the following steps for calculating distance are repeated for each objective
  1. Sort the valid/duplicate adults based on how close they are to that objective using the std::sort function
  2. Add the maximum distance for one objective to the extreme (first and last) adults
  3. Set the normalization value to the worst adult (who has the highest difference between the parameter and the target value)
  4. Add the maximum distance value to all adults meeting the convergence threshold (allowed difference) for this objective
  5. For the rest of the adults not assigned the max distance, it adds the absolute value of both neighbor's objective difference variables divided by the normalization value and subtracted from each other

### mainSort
  - Will sort the inputted adults by [rank-distance](#rankdistancesort) or [rank-rarity](#rankraritysort) depending on which algorithm was specified by the user

### parameterSort
<i>Currently not in use</i>
  - Employs a switch statement to read an inputted objective's goal and sort the inputted vector of Adults by the objective's corresponding sort function

<br>

## planetInfo.cpp/h
  - See the [planet_calculations_info.md](./planet_calculations_info.md) file for detailed information on how planetary body positions are calculated
### PlanetInfo
  - Going backwards in time, this constructor calculates the position of the planet from the endTime to startTime. It fills planetCon with this data starting at the latest point in time to the earliest possible time the spacecraft might have left Earth.
### getCondition/getConditionDev
  - These functions are practically the same - the only difference is getConditionDev can operate on the GPU as well as the CPU.
  - It takes in a time (currentTime) which is the number of seconds BEFORE endTime and will perform an interpolation with the data in planetCon to find the position of planet currentTime seconds before endTime.
### calcIndex
  - Takes in a time value and outputs a corresponding index to be used in the planetCon array
### calc_time
  - Takes in an index and outputs the time corresponded to that index
### getTolData
  - Returns the tolData variable
### getAllPositions
  - Returns the planetCon pointer
### getPlanetPosition
  - Fills the planet variable with planet [elements](#elements)/spacial and velocity information on the planet wanted (either earth or mars)
### interpolate
  - Calculates an approximate [planet element](#elements) for a time between two indexes
### ~PlanetInfo
  - Destructor that removes the dynamic memory allocated to [planetCon](#members-1)
### planetInitial_incremental
  - Determines the next step back from a given planet via [rk4Reverse](#rk4reverse)
### getPlanetSize
  - Returns the size needed for the planetCon array so that [gpuvalues.initialize()](#initialize) can be set up correctly
<br>

## ReferencePoints.cpp/h
- The file contains rank-rarity functions and the [referencePoints class](#referencepoints-class)
### giveRarity
- Calls all the necessary functions to handle rarity assignments for individuals
- Calls [calculateRelCost](#calculaterelcost), then [findAssociatedPoints](#findassociatedpoints), and returns the result of [calculateRarity](#calculaterarity)

### calcNormVector
- Calculates and returns the normal vector for a plane used to calculate the normalization for any number of objectives greater than or equal to two. Called by [calculateRelCost](#calculaterelcost)

### calculateRelCost
- Finds the normalized objective value for each objective of each individual

### findPointDist
- Finds the distance between the passed in point (a vector of doubles) and the normalized adult (normalized in [calculateRelCost](#calculaterelcost))

### findAssociatedPoints
- Finds the closest reference point to each adult in the passed in vector (newAdults)

### calculateRarity
- Assign a unique rarity value to each adult in the allAdults vector. To do this, the algorithm loops until the number of rarities assigned equals the number of adults needed to be assigned a rarity score:
  1. Find the reference point with the least amount of associations
  2. If the reference point does not have an adult associated with it, delete the point from consideration
  3. If the reference point does have at least one adult associated with it, assign the last adult associated with this point the lowest unassigned rarity (lower means more rare) score. 
  4. Remove this adult from consideration

### assignReservedRarity
- If the config defined [reserve rarity](./config_readme.md/#other-variables) score is not 0, then reserve a top rarity score for reserve rarity number of individuals

<br>

## motion_equations.cpp/h
- Contains various functions for [runge_kutta.cpp](#runge_kuttacpph)

### calc_k and calc_kPlanet
- Calculates the corresponding k for the Runge-Kutta computation
- Both functions call [calcRate_r](#calcrate_r), [calcRate_theta](#calcrate_theta), and [calcRate_z](#calcrate_z). calc_k() also calls [calcRate_vr](#calcrate_vr), [calcRate_vtheta](#calcrate_vtheta), and [calcRate_vz](#calcrate_vz) while calc_kPlanet also calls [calcRate_vrPlanet](#calcrate_vrplanet), [calcRate_vthetaPlanet](#calcrate_vthetaplanet), and [calcRate_vzPlanet](#calcrate_vzplanet)

### calcRate_r
- Returns the vr variable from the passed in [elements](#members-1) structure

### calcRate_theta
- Returns vtheta / r. Both variables are from the passed in [elements](#members-1) structure

### calcRate_z
- Returns the vz variable from the passed in [elements](#members-1) structure

### calcRate_vr
- Calculates vrDot

### calcRate_vtheta
- Calculates vthetaDot

### calcRate_vz
- Calculates vzDot

### calcRate_vrPlanet
- Calculates vrPlanetDot

### calcRate_vthetaPlanet
- Calculates vthetaPlanetDot REVIEW

### calcRate_vzPlanet
- Calculates vzPlanetDot

<br>

## output.cpp/h
- Contains functions for displaying genetic algorithm information, and the [output](#output-class) class
### printBestAdults
- Prints the best rank-distance/rarity adult, the best adult for each objective to the terminal, and more generation information. Called in [printGeneration](#printgeneration) and [printFinalGen](#printfinalgen)
### terminalDisplay
- Utility function to display the currently best individual onto the terminal while the algorithm is still running. Called by [printBestAdults](#printbestadults).

<br>

## runge_kutta.cpp/h

### rk4sys
<i>Currently not in use</i>

- Outputs a dynamic array of position and velocity sets. The last entry is the final conditions

###  rk4Reverse
- Calculates the planet's launch date conditions

### rkCalc
- Calculates the k values and gets a new value for y

### rkCalcPlanet
- Calculates the k values

### calc_scalingFactor
- Calculates the scaling factor. Used in [rk4Reverse](#rk4reverse)

### pmLimitCheck
- Looks for an error in step size. Called by [calc_scalingFactor](#calc_scalingfactor)

<br>

## runge_kuttaCUDA.cu/cuh

### callRK
- Fills in the final variables (final position, velocity, etc) for each child. Calls [rk4CUDASim](#rk4cudasim) to simulate the children

### rk4CUDASim
- Simulates the mission for each child

<br>

## calcFourier.cpp/h

### calc_Series

- Calculates the value of the fourier series evaluated at a given time

### calc_gamma

- Calculates the in-plane angle derived from the normalized time and gamma Fourier series. Calls [calc_series](#calc_series)

### calc_tau

- Calculates the in-plane angle derived from the normalized time and tau Fourier series. Calls [calc_series](#calc_series)

### calc_coast

- Calculates whether the spacecraft should accelerate or coast based on the coast Fourier series. Either returns True or False. Calls [calc_series](#calc_series)


</br>

# Classes
## Child
The Child class instances are the individuals generated through crossover that are put through callRK to have their positions, speeds, and error statuses calculated. They should never be sorted, as they lack the metrics needed for sorting (rank and distance/rarity). 
### Child_Members
  - startParams: Holds the information about the starting location of the spacecraft. These are the rkParameters that characterize the spacecraft's trip
  - finalPos: The final position and velocity of the spacecraft when it's at the end of its mission
  - The following variables are used as mission outputs to optimize
    - posDiff: The difference in position between the spacecraft and the center of target at end of the run (in AU)
    - speedDiff: The difference in speed between the spacecraft and the target at the end of the run (in AU/s)
    - horzVelDiff: The horizontal plane angle difference between the final velocities of the spacecraft and the target (in Degrees)
    - vertVelDiff: The vertical plane angle difference between the final velocities of the spacecraft and the target (in Degrees)
    - orbitPosDiff: The difference in position between the spacecraft and the orbital radius of the target at the end of the run (in AU)
    - orbitSpeedDiff: The difference in velocity between the spacecraft and the orbit speed of the target at the end of the run (in AU/s)
    - fuelSpent: The amount of fuel spent by the individual during its simulation (in kg)
    - minMarsDist: the minimum distance the child gets from Mars during the simulation (in AU)
    - orbithChange: The change in angular momentum from an orbital assist (in AU<sup>2</sup>/s)
  - objTargetDiffs: An array that stores the absolute value differences between the child's objective values and the target.
  - normalizedObj: Vector of output normalizations found by comparison to other individuals in the generation
  - progress: How close the individual is to convergence. It is a 0 to 1 scale with 0 being poor and 1 meaning that the individual have completed all objectives
  - avgParentProgress: The average of the individual's two parents' progress values
  - birthday: Keeps track of the generation this individual was created in 
  - stepCount: Counts steps in callRK, needed for orbital missions to keep track of steps
  - The following objectives are used for the child to do multiple simulation cycles (orbital assists)
    - simStartTime: the time from which the simulation will start for the child
    - simStartPos: the position and velocity of the child at the start of the simulation
    - simNum: the number of simulation cycles the child has completed
    - simStatus: the status of the child at the end of its simulation
  - errorStatus: Record of if child is computed correctly, should be set in callRK. If a Child has a NAN in any of its elements or if it gets too close to the Sun or Mars, it will be given an error status to prevent it from passing on its undesirable genes

### Constructors
- Default Constructor: Should never really be used, but set up as a guardrail that sets the functional status to default Child so it could theoretically be detected and removed. Since we are not currently using functional status, perhaps we want to make that an errorStatus instead
- General Child Constructor: The next Child constructor is the only that is normally used. It sets the Child's startParams, using the position of Earth and rkParameters generated randomly or through crossover and mutation. It also sets the status for the Child, along with its parents' progress average and its birthday. It also prepares stepCount to be used elsewhere in the code
- Copy Constructor: This constructor allows you to make a Child that is a clone of another Child. The code might be able to do this without this function being explicitly defined, but this constructor was written as an attempt to fix errors we were getting earlier when we did not have it
### getParameters
- Here, parameter does not refer to rkParameters. Instead, it returns a parameter based on a mission objective
- For example, if POS_DIFF is passed in, this function will return the posDiff of this Child
### getPosDiff
- Calculates and sets the posDiff of a Child, comparing the final position of the Child to that of its target
### getSpeedDiff
- Calculates and sets the speedDiff of a Child, comparing the final speed of a Child to that of its target
### getHorzVelDiff
- Calculates and sets the horzVelDiff of a Child, comparing the final horizontal velocity angle of a Child to that of its target
### getVertVelDiff
- Calculates and sets the vertVelDiff of a Child, comparing the final vertical velocity angle of a Child to that of its target
### getOrbitPosDiff
- Calculates and sets the orbitPosDiff of a Child, comparing the final position of the Child to that of its target and the orbital radius the mission is trying to achieve
### getOrbitSpeedDiff
- Calculates and sets the orbitSpeedDiff of a Child, comparing the final speed of the Child to the final speed of its target and factoring in the fact that we are trying to achieve a certain orbital speed
### getProgress
- Calculates and sets the progress of an individual
- Process progress is calculated by:
  1. For each objective, find the upper and lower objective bounds of allowed difference relative to the target value and determine which one is closer to the individual's parameter. If the objective is converged, add 1 to the calculated progress (originally set to 0) and if it's not converged, add 1 plus the absolute value of the difference between the closest bound and the objective value, divided by the closest bound
  2. Once the algorithm has gone through each objective, it divides the number of mission objectives by this calculated progress to find that individual's progress
- If a Child is neither VALID nor a duplicate, its progress is set to 0 because this individual is not useable

## Adult 
The Adult class is built upon the [Child class](#child) and extends it. Individuals that are Adults should have already been through callRK as a Child, and should never be passed through this function again. An Adult inherits all the properties of a Child, but also contains rank and distance/rarity, which allows them to be sorted. Unlike a Child, which should never be sorted, Adults are made to be sorted and they are the individuals that are used to create the next generation of Children. Note: Info on this class does not include some of the unused parameter sort functions
### Adult_Members
- rank: rank is a metric that tells you how an individual compares to other individuals in terms of how well it meets the objectives. Individuals that have a lower rank are better than individuals with a higher rank. For more information about how rank is calculated, please refer to the section about [giveRank](#giverank)
- distance: distance is a metric that tells you how different from other individuals an Adult is. It is a measure of genetic diversity and higher distances (more unique individuals) are favored in sorting. For more information on how an Adult gets a distance, refer to the [giveDistance](#givedistance) section
- rarity: rarity is similar to distance in that is measured how unique it is compared to other individuals. Lower rarity scores mean the adult is considered more unique. For more information about how rarity is calculated, please refer to the section about [calcRarity](#calculaterarity)
- associatedRefPoint: Will hold the index of the reference point that is closest to this Adult when using the rank-rarity algorithm  
### Constructors
- Default Constructor: Sets rank to INT_MAX, distance to -1, rarity to INT_MAX, and associatedRefPoint to 0
- General Adult Constructor: Takes in a child and creates an adult from it
### bestProgress
- This function is passed into std::sort which sorts by progress.
### worstProgress
- This function is passed into std::sort which sorts by worst progress.
### lowAssocPoint
- This function is passed into std::sort which sorts by which reference point an individual is associate with.
### dominationCheck
- In essence, an individual "dominates" another individual if and only if it is better at one objective and at least as good as the other individual in all other objectives
- If both Adults have an objective within the goal difference (domination threshold), this parameter will be ignored in domination
  - If all the parameters are ignored, neither individual dominates the other individual
  - If all the parameters are approximately equal except one, then 1 individual will dominate the other
- If an Adult has ANY parameter where it performs better than the other individual, it will not be dominated by that individual, even if the second individual outperforms it on every other parameter
  - Example: Minimization of posDiff, speedDiff, & tripTime (NOTE: these are not real values). Adult A has posDiff: 1.0e-9; speedDiff: 3.5e-7; tripTime: 47000000 and Adult B has posDiff: 1.0e-8; speedDiff: 1e-10; tripTime:55000000. Adult A's posDiff and tripTime are both better than Adult B's, but Adult B's speedDiff is better than Adult A's, so they are co-dominant. Neither Adult dominates the other. 
  - For a hand-calculated example of this code, refer to the Calculations for testing_genetics PDF in the 2022-Summer file on Teams.
- Domination is used to determine rank. 
### rankSort
- This function is passed into std::sort which sorts by rank, with the lowest rank being favored over higher ranks.
### rankDistanceSort
- This function is passed into std::sort which sorts first by rank, with the lowest rank being favored over higher ranks. Then, if two individuals are in the same rank, the one with the greatest distance (the most unique individual) is favored.
### rankRaritySort
- This function is passed into std::sort which sorts first by rank, with the lowest rank being favored over higher ranks. Then, if two individuals are in the same rank, the one with the lowest rarity (the most unique individual) is favored.
### duplicateCheck
- Compares two adults and returns True is a duplicate is found. A duplicate is an individual who has all the same parameters (within a tolerance) as another individual.
### findDuplicates
- Takes in two Adult vectors, looks for duplicates, and sets their status accordingly


## GPUMem
### Members
  - absTol: Used to be set to rk_tol
  - numThreads: Used to hold num_individuals
  - *devGeneration: Device pointer to allocate space for generation
  - *devAbsTol: Device pointer to allocate space for tolerance
  - *devCConstant: Device pointer to allocate space for cConstants
  - *devMarsLaunchCon: Device pointer to allocate space for marsLaunchCon, the array of positions of Mars over time
  - *devTime_steps: Device pointer to allocate space for each time step
  - *devY_steps: Device pointer to allocate space for the position and velocity of each time step
  - *devGamma_steps: Device pointer to allocate space for the gamma value of each time step
  - *devTau_steps: Device pointer to allocate space for the tau value of each time step
  - *devAccel_steps: Device pointer to allocate space for the acceleration value of each time step
  - *devFuel_steps: Device pointer to allocate space for the fuel value of each time step
### Initialize
- Initializes the device parameters and allocates the space on the GPU
### Free
- Deallocates memory stored on the GPU
## PlanetInfo(Class)
- See the [planet_calculations_info.md](../Cuda/Planet_calculations/planet_calculations_info.md), the section of this document on [PlanetInfo](#planetinfocpph) or read the comment headers on the planetInfo.h file
## Elements
### Members
- These are the 6 main "elements" of the spacecraft: the 3 positional cylindrical coordinates: r, theta, and z; and the spacecraft's velocity in those 3 directions: vr, vtheta, and vz. These 3 components are initialized in elements.cpp and elements.h
### Constructors
- Default Constructor: Sets all members to 0
- General Constructor: Takes in all components as arguments and sets the member variables accordingly
### operator+
- Adds the like member variables from two element structs
### operator-
- Subtracts the like member variables from two element structs
### operator*
- Takes in a scalar value and multiplies each member variable by this value
### operator/
- Takes in a scalar value and divides each member variable by this value
### compare
- Returns a percentage signifying how different the two elements are
### operator<<
- Returns all 6 member varaibles. Allows for writing to a file. (or printing to a terminal)

## rkParameters
### Description
- These are the 19 "genes" of the spacecraft. This includes the 6 parameters in elements, triptime, beta, zeta, alpha, and the coefficients for gamma (7), tau (3), and coast (5). These genes are initialized in rkParameters.cpp and rkParameters.h.   REVIEW HOW MANY GENES??
### Members
- y0: An [elements struct](#elements) to hold initial position and velocity elements
- coeff: A [coefficients struct](#coefficients) to hold the arrays for the fourier series
- tripTime: final time of the simulation
- alpha: position portion of the launch angles when leaving earth's SOI CORRECT???
- beta: in plane velocity portion of the launch angles when leaving earth's SOI
- zeta: out-of-plane velocity portion of the launch angles when leaving earth's SOI
### Constructors
- Default Constructor: Sets all members to 0 (except coeff since it was giving us errors when setting it to NULL)
- General Constructors: These constructors take in various combinations of arguements to set some (or all) of the member variables
### randomParameters
- Returns a randomized rkParameters object
### compare
- Returns a percentage signifying how different the two rkParameter objects are
### parametersRK4Simple
<i>Currently not in use</i>
- Calls rk4Simple which doesn't seem to exist

## objective
### Description
- Holds information on each objective. This file also holds the parameterGoals enumeration
### Members
- name: The name the program will use for the objective
- goal: the parameter goal (See description) of the objective
- target: the target/convergence threashold of the objective (defined in the [mission config readme](./mission_config_readme.md))
- allowedDifference: the allowed difference for the objective (defined in the [mission config readme](./mission_config_readme.md))
- goalDifference: the goal difference for the objective (defined in the [mission config readme](./mission_config_readme.md))
- equateTolerance: the equate tolerance for the objective (defined in the [mission config readme](./mission_config_readme.md))
### Constructors
- Default Constructor: Sets the goal as UNKNOWN (based on the parametersGoals enumeration mentioned in the description)
- General Constructor: Takes in arguements for each member variable and sets the member variables accordingly

## ReferencePoints (class)
### Description
- Handles the creation of reference points
### Members
- points: A 2D vector which will store all reference points and their individual components (values for each objective so that each point exists at a different location in objective space)
- objWorst: will store the worst individual when [calculating the relative cost](#calculaterelcost) (normalization of each objective)
- objBest: will store the best individual when [calculating the relative cost](#calculaterelcost) (normalization of each objective)
### Constructors
- Default Constructor: Should not be used and is not defined
- General Constructor: Will calculate all of the reference points based on the number of objectives and the [divisions per objective](./config_readme.md/#other-variables) by calling [addPoint](#addpoint). When this constructor is finished, we should have evenly distributed reference points across each normalized objective space.
### addPoint
- Will add a new point to the [points](#description-2) vector 

## output (class)
### Members
- baseFolder: Will store the overall folder where each individual run's files are stored
- outputPath: Will store the folder address for the individual run's output files
### Constructor
- Default (and specific) Constructor: sets the base folder and calls [prepareFolders](#preparefolders) and [initializeGenPerformance](#initializegenperformance) to set data collection up
### printGeneration
- Handles the printing and recording of data during a run
### printFinalGen
- Once the run has finished, record the final generation's information, including the needed information to plot the best spacecraft's path.
### initializeGenPerformance
- Create the genPerformance .csv file for later output by [recordGenPerformance](#recordgenerationperformance)
### initializeSimpleGenPerformance
<i>Currently not in use</i>

- Creates a simpler version of genPerformance for later output

### prepareFolders
- Creates the current run's folder for data output

### recordGenerationPerformance
- Appends the current generation's data to the genPerformance .csv file

### recordGenSimple
<i>Currently not in use</i>

- Appends simple generational data to the simpleGenPerformance .csv file

### recordAllIndividuals
- records all individuals of the current generation, along with their data.

### reportRun
- Prints general information about a run when it is over

### trajectoryPrint
- gets trajectory data for the finished run and puts it in the orbitalMotion and finalOptimization .bin files

### errorCheck
- Records error in energy conservation due to thrust calculations

### recordEarthData
- Creates a earthCheckValues .csv file to hold rows of [element](#elements) info on earth with a time stamp on each row

### recordMarsData
- Creates a marsCheckValues .csv file to hold rows of [element](#elements) info on mars with a time stamp on each row

### recordReferencePoints
- Print a run's used reference points

## coefficients

### Members
- gammaSize: Sets the size of the gamma array
- gamma: A fourth order Fourier series array
- tauSize: Sets the size of the tau array
- tau: A first order Fourier series array
- coastSize: Sets the size of the coast array
- coast: A second order Fourier series array
- coastThreshold: If the value is above the threshold, acceleration occurs. When the value is below, coasting occurs

### initCoefficient
<i> Doesn't seem to be used </i>

- Initializes coefficients?

### operator<<
<i> Doesn't seem to be used </i>

- Returns the values in the gamma, tau, and coasting member variables

## thruster

### Members
- P0: Initial power in
- type: The type of thruster to use. Uses a enumeration called THRUST_TYPE defined in the same file
- coastThreshold: Config defined coast threashold determines when the spaceship should coast and when it should accelerate
- Maximum power and flow rates for different variables

### Constructor
- Default Constructor: Uses [cudaConstants](#configcpph) to determine thruster type

### calc_eff
- Evaluates the spacecraft's efficiency for a certain iteration

### calc_m_Dot
- Evaluates the spacecraft's fuel flow rate (mDot) for a certian iteration

### calc_accel
- Evaluates the acceleration of the spacecraft