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
- [Classes](#classes)


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

<br>

## Anneal.cpp/h
### changeAnneal
  - Calculates the anneal for the next generation
    - Anneal affects the strength of mutations to children parameters, which allows for smaller changes to input parameters when close to solving the objectives
  - Current implimentation relies on a step method 
    - Each step sets the current anneal to a multiplier of a intial anneal
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
  - This differs from crossOver_wholeRandom because the parent for the gamma, tau, and coast parameters to one parent

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
    1. fillParents; this selects the parents to be used when generating children
    2. makeChildren; this creates the new set of children
    3. callRK; this passes the children along to be simulated
    4. convertToAdults; this fills newAdults with the simulated children

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
  2. Each individual is assigned a rank depending on how many times they were dominated
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
  - 

<br>
</br>

# Classes
