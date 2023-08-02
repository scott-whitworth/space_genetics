<h1>Code Overview (CS)</h1>
<i>This document is intended to explain the basic purpose behind the functions and classes used throughout the code, specifically though a CS perspective</i>
<br>
<i>Last Updated: Aug 2022</i>
<br>
</br>

# Table of Contents
- [Table of Contents](#table-of-contents)
- [Surface Level Info](#surface-level-info)
- [Files](#files)
  - [Optimization.cu](#optimizationcu)
    - [main](#main)
    - [optimize](#optimize)
    - [checkTolerance](#checktolerance)
    - [calculateGenerationValues](#calculategenerationvalues)
  - [config.cpp/h](#configcpph)
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
- [Classes](#classes)
  - [Child](#child)
    - [Child_Members](#child_members)
    - [Constructors](#constructors)
    - [Unit_Test_Constructors](#unit_test_constructors)
    - [getParameters](#getparameters)
    - [getPosDiff](#getposdiff)
    - [getSpeedDiff](#getspeeddiff)
    - [getOrbitPosDiff](#getorbitposdiff)
    - [getOrbitSpeedDiff](#getorbitspeeddiff)
    - [getProgress](#getprogress)
  - [Adult](#adult)
    - [Adult_Members](#adult_members)
    - [dominationCheck](#dominationcheck)
    - [rankSort](#ranksort)
    - [rankDistanceSort](#rankdistancesort)
  - [GPUMem](#gpumem)
    - [Members](#members)
    - [Initialize](#initialize)
    - [Free](#free)
  - [PlanetInfo(Class)](#planetinfoclass)
  - [Elements](#elements)
    - [Description](#description)
  - [rkParameters](#rkparameters)
    - [Description](#description-1)



<br>
</br>

# Surface Level Info

<br>
</br>

# Files

## Optimization.cu
### main
  - Function which starts the algorithm
  - Gathers CUDA device properties, reads the inputted config variables, calculates the possible positions for the Earth and Mars, and sets the run's base rng seed 
  - Runs a for loop which will call the optimize function for the number of runs specified in the config files
    - Recalculates the potential positions for planetary bodies each time through the loop

### optimize
  - This function facilitates finding the optimal path for the input mission
  - It initially creates the 0th generation, made up of adults that are either imported or randomly generated
  - It then starts a do while loop which ends when either the generation limit has been reached or it has found an optimized answer
    - One run through this loop consitutes one generation
  - First step in the loop is to generate new childeren from existing adults
  - Then, each individual is assigned a rank and distance determined by being compared to all other individuals
  - Then, generation-wide values are calculated for this generation
    - These values are primarily min/max/avg of outputs from the individuals' simulations
  - It will then calculate the anneal for the next generation
  - Output files are then created or appended to with up-to-date information
  - Finally, there is a check to see if the best individual(s) have found an answer
  - If an answer has been found, final outputs are created

### checkTolerance
  - Checks the top x number of individuals to see if they have found an optimized answer
    - x is defined in config files
  - It dynamically checks parameters depending on the number of input objectives
  - It returns true if the top x individuals have found an answer

###  calculateGenerationValues
  - Calculates generation-wide values which are passed to functions which write output files
  - List of values calculated:
    - Average, minimum, and maximum distance values
    - Average and oldest ages and birthdays 
    - Average and best values for each objective

<br>

## config.cpp/h
- Sets us the constants that will be used throughout the code
- See [config_readme.md](Config_Constants/config_readme.md) for more information on this section and for information on the config files themselves
<br>

## Anneal.cpp/h
### changeAnneal
  - Calculates the anneal for the next generation
    - Anneal affects the strength of mutations to children parameters, which allows for smaller changes to input parameters when close to solving the objectives
  - Current implementation relies on a step method 
    - Each step sets the current anneal to a multiplier of a initial anneal
    - The step is chosen based on the progress of the best individual
    - The higher the step, the lower the anneal is set to
    - When the progress of the best individual gets high enough, the anneal calculation method switches to a method which employs exponential decay
  - The base current anneal set by the steps is modified by a sinusoidal-based multiplier
    - As of writing, the current anneal is modified by up to +/- 70% the base value of the current anneal
      - the modifier value is determined by a sinusoidal curve based on the current generation
  - Finally, it will make sure the new anneal never dips below the minimum anneal

<br>

## ga_crossover.cpp/h
### maskValue enum
  - Used to generate masks for generating children
  - Determines if a child's input parameter comes from parent 1 (1), parent 2 (2), the average of the parents (3), or a weighted average of the parents (4)
  - An array with the size of the number of input parameters is filled with instances of this enum; this is used to determine the origin of each of a child's input parameters

### crossOver_randHalf
  - Sets the the mask's indexes to be either parent 1 or parent 2
  - In this function, a cross index is randomly generated
  - Any index before the cross index is assigned to parent 1 and any index after is assigned to parent 2

### crossOver_oneParent
  - Assigns all indexes of the mask array to parent 1

### crossOver_wholeRandom
  - Randomly assigns all indexes of the mask array to either parent 1 or parent 2

### crossOver_bundleVars
  - Randomly assigns the indexes of the mask array to either parent 1 or parent 2
  - This differs from crossOver_wholeRandom because the parent for all the gamma, tau, or coast parameters are just one parent
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
  - Takes in two crossover mask arrays and copies all of the enum values in the first/maskIn array to the second/maskOut array
  - <i>Note: currently unused</i>

### getRand
  - Generates a random number between the positive and negative value of an inputted maximum value
    - For example, if the argued maximum is 10, it will return a random value between 10 and -10
  - Used for generating input parameter mutations

### mutateMask
  - Fills in a mutation mask for a new child
  - Mask, with the size of the number of input parameters, is marked true for all input parameters to mutate
  - The function generates a random number between 0 and 1 and will mutate a parameter if the number is larger than the mutation rate defined in genetic.config
    - If the number is smaller, the function will stop marking parameters for mutation
  - If the function determines a parameter needs to be mutated, it will randomly select a parameter to mutate
    - It will keep randomly choosing parameters to mutate until it finds one that hasn't been marked for mutation
  - Function used only by the mutate() function

### mutate
  - Mutates children newly created from crossover
  - Uses mutateMask to determine which of the child's parameters to mutate
  - To mutate a parameter it uses getRand(), the current anneal, and the input parameter mutate scales from genetic.config to modify input parameters marked for mutation by a random +/- value
    - The anneal multiplier means the mutations will get smaller as the best individuals get closer to solving the objectives

### generateNewChild
  - Generates a set of input parameters for a new child based on a crossOver mask
  - The parameters come from two adults passed in
  - The mask is an array with indexes corresponding to each input parameter. Each index either has enum values of Partner1, Partner2, Avg, or Avg_Ratio, which will determine the origin of each of the child's input parameters
  - If the mask index for a parameter has Parent1 as the value, the parameter will be pulled from the first adult passed in
  - If the mask index for a parameter has Parent2 as the value, the parameter will be pulled from the second adult passed in
  - If the mask index for a parameter has Avg as the value, the parameter will be set as the 50/50 average of the parents' parameters
  - Finally, if the mask index for a parameter has Avg_ratio as the value, the parameter will be set as a weighted average of the parents' parameters
    - The weight of the average is randomly generated at the beginning of the function
  - After a set of parameters is generated, it will be mutated using the mutate function

### generateChildrenPair 
  - Creates two new children using generateNewChild and appends them to an array with new children
  - Flips the crossOver children after generating the first child, so two unique children should be generated
    - Exceptions may be if the parents are duplicates or if the average_ratio randomly generates a 50/50 average
  - Function will exit after generating the first child if the limit for children generation has been reached

### generateChildrenFromCrossover
  - Generates a number of children by the genetic crossover method
    - The number of children generated comes from an inputted value 
  - Creates parent pairs by iterating through a shuffled list of the potential parents 
    - Note: The potential parents themselves are not shuffled, rather a vector with indexes corresponding to each future parent is shuffled
  - As of writing, it generates 6 children from three crossover methods
    - The used methods are crossOver_wholeRandom, crossOver_average, and crossOver_bundleVars
    - The same two parents generate all 6 new children
    - The number of children a pair of parents generate is flexible
  - If all parents have created children and children still need to be generated, the potential parents will be reshuffled, allowing for new parent pairs

<br>

## geneticAlgorithm.cpp/h
### newGeneration
  - Calls all functions necessary for generating a new set of children which will fill the newAdults vector
  - Functions called are:
    1. fillParents- this selects the parents to be used when generating children
    2. makeChildren- this creates the new set of children
    3. callRK- this passes the children along to be simulated
    4. convertToAdults- this fills newAdults with the simulated children

### fillParents
  - Selects the top adults from the previous generations to be parents
    - The number of parents considered comes from the survivor_count variables in genetic.config
  - If it is generation 0, parents are selected from the first adults who were randomly generated
  - Potential parents will only include adults with error statues of valid or duplicate

### makeChildren
  - Currently, it only calls the generatingChildrenFromCrossover function, making sure the number of children generated equals num_individuals, which is defined in genetic.config

### convertToAdults
  - This first completes final output calculations and error checks for the new children
    - If it detects a child has nans for final position or velocity, it will assign an error to the child
    - Otherwise, it will calculate outputs which need their own dedicated function
  - Then, all new children are pushed to the newAdults vector (which is emptied at the beginning of the function)

### createFirstGeneration
  - This function creates a number of new children which are either randomly generated or pulled from a document and passes them to the firstGeneration function
  - Note: as of writing, the functionality of pulling initial children from a document has not been tested since 2020 and may not work

### firstGeneration
  - Passes the children generated from createFirstGeneration to the simulation using callRK
  - Converts the simulated initial children into adults and pushes them into the oldAdults vector

### preparePotentialParents
  - Function facilitates the calculation of rank and distance for all individuals
    - This allows for the right parents to be chosen next generation and for the best individual to be printed and checked for convergence
  - It will first determine which individuals are duplicates using findDuplicates
  - Then it will add all valid adults to the allAdults vector using addToAllAdults
    - Note: the change which will get rid of or keep duplicates is within this function
  - Then it will calculate the rank and distance for all individuals in allAdults using the respective functions (giveRank and giveDistance)
  - Finally, it will replace the individuals in the oldAdults vector with the top half rank-distance individuals from allAdults

### addToAllAdults
  - Goes through all individuals in the inputted adult vector and checks their error statuses
    - If an adult's error status is not valid or duplicate, one is added to an error tracking variable used for reporting
    - If an adult's error status is duplicate, one is added to a duplicate tracker used for reporting and the adult is appended to allAdults
      - Note: if you want to remove duplicates when they are detected, comment out the line where duplicate adults are pushed to allAdults
    - If an adult's error status is valid, it is pushed to allAdults

<br>

## Sort.cpp/h
### giveRank
  - Assigns a rank to all valid (or duplicate) individuals
  - A rank can be thought of as more of tiers, where multiple individuals can occupy the same rank
    - However, when solving for only one parameter, rank works as expected, with each unique individual occupying a unique rank
  - Assigning a rank occurs in two stages:
  1. Each individual is compared against every other individual passed in and there is a check to see if one individual dominates the other
      - See the section on the adult class to learn about how an individual dominates another
      - For each individual, it is tracked how many times they were dominated and what other individuals they have dominated
      - At the end of the comparisons, any individual who was never dominated is added to a "front" vector
  1. Each individual is assigned a rank depending on how many times they were dominated
       - This process utilizes a while loop which runs until all individuals have a rank; the loop is:
         1. Go through all individuals in the front vector and pull up the individuals they dominated
         2. For each of the dominated individuals, subtract one from their counter for how many times they were dominated
         3. If an individual's domination counter reaches 0, add them to the "newFront" vector and assign them their rank (rank = times though this loop + 1)
         4. Once the front vector has been iterated through, replace the front vector with the newFront vector and clear the newFront vector
            - If newFront is empty during this step, end the loop instead of repeating

### giveDistance
 - Assigns a distance to all valid (or duplicate) individuals
 - Distance measures how diverse an individual is compared to its neighbors
   - It is calculated by dividing the difference between the parameters of the individuals immediately better and worse than the focus individual by a normalization value
     - The normalization value is always the largest value for the parameter within the population of individuals
 - First counts the number of valid or duplicate individuals so only they can be compared to later
 - Makes sure all distances are reset to 0 before calculating new distances
 - Note: the following steps for calculating distance are repeated for each objective
  1. Sort the valid/duplicate adults in the correct order for the objective using the parameterSort function
  2. Adds the maximum distance for one objective to the extreme (first and last) adults
  3. Checks if the goal is a minimum or maximum
    - If it is a minimization, it sets a normalization value to the last individual's variable
    - If it is a maximization, it sets the normalization value to the variable of the first individual who hasn't met the convergence tolerance
  4. It adds the maximum distance value to all adults meeting the convergence threshold for this objective
  5. For the rest of the adults not assigned the max distance, it adds the value of the difference between the adult's neighbors divided by the normalization value to the adult's distance

### parameterSort
  - Employs a switch statement to read an inputted objective's goal and sort the inputted vector of Adults by the objective's corresponding sort function

<br>

## planetInfo.cpp/h
  - See the [planet_calculations_info.md](../Cuda/Planet_calculations/planet_calculations_info.md) file for detailed information on how planetary body positions are calculated
### PlanetInfo
  - Going backwards in time, this function calculates the position of the planet from the endTime to startTime. It fills planetCon with this data starting at the latest point in time to the earliest possible time the spacecraft might have left Earth.
### getCondition/getConditionDev
  - These functions are practically the same - the only difference is getConditionDev can operate on the GPU as well as the CPU.
  - It takes in a time (currentTime) which is the number of seconds BEFORE endTime and will perform an interpolation with the data in planetCon to find the position of planet currentTime seconds before endTime.

<br>
</br>

# Classes
## Child
The Child class are the individuals generated through crossover that are put through callRK to have their positions, speeds, and error statuses calculated. They should never be sorted, as they lack the metrics needed for sorting (rank and distance). 
### Child_Members
  - startParams: Holds the information about the starting location of the spacecraft. These are the rkParameters that characterize the spacecraft's trip
  - finalPos: The final position and velocity of the spacecraft when at the end of its mission
  - posDiff: The difference in position between spacecraft and center of target at end of run (in AU)
  - speedDiff: The difference in speed between spacecraft and target at end of run (in AU/s)
  - orbitPosDiff: The difference in position between the spacecraft and the orbital radius of the target at the end of the run (in AU)
  - orbitSpeedDiff: The difference in velocity between the spacecraft and the orbit speed of the target at the end of the run (in AU/s)
  - fuelSpent: The amount of fuel spent by the individual during its simulation (in kg)
  - progress: How close the individual is to convergence. It is a 0 to 1 scale with 0 being poor and 1 meaning that the individual have completed all objectives
  - avgParentProgress: The average of the individual's two parents' progress values
  - birthday: Keeps track of the generation this individual was created in 
  - stepCount: Counts steps in callRK, needed for orbital missions to keep track of steps
  - errorStatus: Record of if child is computed correctly, should be set in callRK. If a Child has a NAN in any of its elements or if it gets too close to the Sun or Mars, it will be given an error status to prevent it from passing on its undesirable genes

### Constructors
- Default constructor: Should never really be used, but set up as a guardrail that sets the functional status to default Child so it could theoretically be detected and removed. Since we are not currently using functional status, perhaps we want to make that an errorStatus instead
- General Child Constructor: The next Child constructor is the only that is normally used. It sets the Child's startParams, using the position of Earth and rkParameters generated randomly or through crossover and mutation. It also sets the status for the Child, along with its parents' progress average and its birthday. It also prepares stepCount to be used elsewhere in the code
- Copy Constructor: This constructor allows you to make a Child that is a clone of another Child. The code might be able to do this without this function being explicitly defined, but this constructor was written as an attempt to fix errors we were getting earlier when we did not have it
### Unit_Test_Constructors
- A variety of constructors were created specifically for unit tests. They are within a #ifdef, so they should not be accessible in the normal code, as they could potentially cause major errors if they were used anywhere aside from unit testing. For the specifics of what these constructors are used for, refer to [child.h](Genetic_Algorithm/child.h).
### getParameters
- Here, parameter does not refer to rkParameters. Instead, it returns a parameter based on a mission objective.
- For example, if MIN_POS_DIFF is passed in, this function will return the posDiff of this Child.
- This function helps sort individuals by different objectives, by pulling the part of Child that corresponds to a mission objective.
### getPosDiff
- Calculates and sets the posDiff of a Child, comparing the final position of the Child to that of its target.
### getSpeedDiff
- Calculates and sets the speedDiff of a Child, comparing the final speed of a Child to that of its target.
### getOrbitPosDiff
- Calculates and sets the orbitPosDiff of a Child, comparing the final position of the Child to that of its target and the orbital radius which the mission is trying to achieve.
### getOrbitSpeedDiff
- Calculates and sets the orbitSpeedDiff of a Child, comparing the final speed of the Child to the final speed of its target and factoring in the fact that we are trying to achieve a certain orbital speed.
### getProgress
- Calculates and sets the progress of an individual. 
- Progress is calculated by finding the a ratio between the goal/convergence threshold and the actual value an individual possesses for a certain objective. If the individual has met an objective, a 1 is added for that metric, otherwise it will be a ratio of the threshold and actual values for an objective. Once the actual vs. goal ratio has been calculated for every objective, the progress should be a number greater than or equal to the number of objectives. Then the number of objectives is divided through by the number calculated for progress so it is a number between 0 and 1, where 1 means an individual has converged.
  - The progress calculation is done a little differently depending on whether the goal is a minimization or maximization. If the goal is a minimization, the actual value is divided by the threshold. If the goal is a maximization, the threshold is divided by the actual value.
- If a Child is neither VALID nor a duplicate, its progress is set to 0 because this individual is not useable.

## Adult 
The class Adult is built upon the Child class and extends it. Individuals that are Adults should have already been through callRK as a Child, and should never be passed through this function again. An Adult inherits all the properties of a Child, but also contain rank and distance, which allows them to be sorted. Unlike a Child, which should never be sorted, Adults are made to be sorted and they are the individuals that are used to create the next generation of Children.
### Adult_Members
- rank: rank is a metric that tells you how an individual compares to other individuals in terms of how well it meets the objectives. Individuals that have a lower rank are better than individuals with a higher rank. For more information about how rank is calculated, please refer to the section about [giveRank](#giverank). 
- distance: distance is a metric that tells you how different from other individuals an Adult is. It is a measure of genetic diversity and higher distances (more unique individuals) are favored in sorting. For more information on how an Adult gets a distance, refer to the [giveDistance](#givedistance) section.
### dominationCheck
- In essence, an individual "dominates" another individual if and only if it is better (or at least as good as) another individual for EVERY objective.
- If two Adults have an objective within equate tolerance or below the domination threshold, this parameter will be ignored
  - If all the parameters are ignored, neither individual dominates the other individual
  - If all the parameters are approximately equal except one, then 1 individual will dominate the other
- If an Adult has ANY parameter where it performs better than the other individual, it will not be dominated by that individual, even if the second individual outperforms it on every other parameter
  - Example: Minimization of posDiff, speedDiff, & tripTime (NOTE: these are not real values). Adult A has posDiff: 1.0e-9; speedDiff: 3.5e-7; tripTime: 47000000 and Adult B has posDiff: 1.0e-8; speedDiff: 1e-10; tripTime:55000000. Adult A's posDiff and tripTime are both better than Adult B's, but Adult B's speedDiff is better than Adult A's, so they are co-dominant. Neither Adult dominates the other. 
  - For a hand-calculated example of this code, refer to the Calculations for testing_genetics PDF in the 2022-Summer file on Teams.
- Domination is used to determine rank. 
### rankSort
- Individuals that are either VALID or DUPLICATES are sorted by their rank, with the lowest rank being favored over higher ranks.
### rankDistanceSort
- Individuals that are either VALID or DUPLICATES are sorted by their rank, with the lowest rank being favored over higher ranks. If two individuals are in the same rank, the one with the greatest distance (the most unique individual) is favored.

## GPUMem
### Members
  - stepSize: Not currently in use 
  - absTol: Used to be set to rk_tol
  - numThreads: Used to hol num_individuals
  - *devGeneration: Device pointer to allocate space for generation
  - *devTimeInitial: Device pointer to allocate space for timeInitial (set to zero)
  - *devStepSize: Device pointer to allocate space for the step size (not currently in use)
  - *devAbsTol: Device pointer to allocate space for tolerance
  - *devCConstant: Device pointer to allocate space for cConstants
  - *devMarsLaunchCon: Device pointer to allocate space for marsLaunchCon, the array of positions of Mars over time
### Initialize
- Initializes the device parameters and allocates the space on the GPU
### Free
- Deallocates memory stored to GPU
## PlanetInfo(Class)
- See the [planet_calculations_info.md](../Cuda/Planet_calculations/planet_calculations_info.md), the section of this document on [PlanetInfo](#planetinfocpph) or read the comment headers on the planetInfo.h file
## Elements
### Description
- These are the 6 main "elements" of the spacecraft being the 3 positional cylindrical coordinates: r, theta, and z; and the spacecraft's velocity in those 3 directions: vr, vtheta, and vz. These 3 components are initialized in elements.cpp and elements.h and used as double throughout the code. They are calculated within child after callRK.
## rkParameters
### Description
- These are the 25 (26?) "genes" of the spacecraft. This includes the 6 parameters in elements, triptime, beta, zeta, alpha, and the coefficients for gamma (7), tau (3), and coast (5). These genes are initialized in rkParameters.cpp and rkParameters.h. Other than elements, these are all calculated when callRK is called, and are mutated in ga_crossover.cpp. 