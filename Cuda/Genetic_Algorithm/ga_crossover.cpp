#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>

#include "..\Genetic_Algorithm\adults.h"

// Global enumeration for the different mask values instead of 1,2,3 for better readibility and clarity of value meaning
enum maskValue {
    PARTNER1 = 1,
    PARTNER2,
    AVG,
};

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Creates a random bifurcation mask, currently not in use
void crossOver_randHalf(std::vector<int> & mask, std::mt19937_64 & rng) {
    int crossIndex = rng() % (OPTIM_VARS-1);
    //cout << "Random Index: " << crossIndex << endl;
    for (int i = 0; i < OPTIM_VARS; i++) {
        if (i > crossIndex) {
            mask[i] = PARTNER2;
        }
        else {
            mask[i] = PARTNER1;
        }        
    }
    return;
}

// Sets the entire mask to be PARTNER1 for length OPTIM_VARS, 
// allows a crossover where no mixing occurs
// currently not in use
void crossOver_oneParent(std::vector<int> & mask) {
    for (int i = 0; i < OPTIM_VARS; i++) {
        mask[i] = PARTNER1;
    }
    return;
}

// Creates a random mask
void crossOver_wholeRandom(std::vector<int> & mask, std::mt19937_64 & rng) {
    for (int i = 0; i < OPTIM_VARS; i++ ) {
        if (rng() % 2) { //Coin flip, either 1/0
            mask[i] = PARTNER2;
        }
        else {
            mask[i] = PARTNER1;
        }
    }
    return;
}

//Generates crossover mask that maintains paramter relationships
void crossOver_bundleVars(std::vector<int> & mask, std::mt19937_64 & rng) {
    int p_gamma = 1 + rng() % 2; // partners for each variable are randomly chosen between 1 and 2
    int p_tau = 1 + rng() % 2;
    int p_coast = 1 + rng() % 2;
    int p_triptime = 1 + rng() % 2;
    int p_alpha = 1 + rng() % 2;
    int p_beta = 1 + rng() % 2;
    int p_zeta = 1 + rng() % 2;

    // Copying random assignment to all of paramters
    // Gamma values
    for (int i = GAMMA_OFFSET; i < (GAMMA_OFFSET + GAMMA_ARRAY_SIZE); i++) {
        mask[i] = p_gamma;
    }
    // Tau values
    for (int i = TAU_OFFSET; i < (TAU_OFFSET + TAU_ARRAY_SIZE); i++) {
        mask[i] = p_tau;
    }
    // Coast values
    for (int i = COAST_OFFSET; i < (COAST_OFFSET + COAST_ARRAY_SIZE); i++) {
        mask[i] = p_coast;
    }

    mask[TRIPTIME_OFFSET] = p_triptime;
    mask[ALPHA_OFFSET] = p_alpha;
    mask[BETA_OFFSET] = p_beta;
    mask[ZETA_OFFSET] = p_zeta;

    return;
}

// Sets the entire mask to be AVG for length OPTIM_VARS
void crossOver_average(std::vector<int> & mask) {
    for (int i = 0; i < OPTIM_VARS; i++) {
        mask[i] = AVG;
    }
    return;
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Utility to flip the polarity of a mask
void flipMask(std::vector<int> & mask) {
    for (int i = 0; i < OPTIM_VARS; i++) {
        if (mask[i] == PARTNER1) {
            mask[i] = PARTNER2;
        }
        else if (mask[i] == PARTNER2) {
            mask[i] = PARTNER1;
        }
        // If mask[i] is neither partner1 nor partner2 (must be avg then), leave it be
        // To handle the average enumeration
    }
    return;
}

// Copy contents of maskIn into maskOut of size OPTIM_VARS
void copyMask(int * maskIn, int * maskOut) {
    for (int i = 0; i < OPTIM_VARS; i++) {
        maskOut[i] = maskIn[i];
    }
}

// Utility function for mutate() to get a random double with high resolution
// The range of the random number is from -max to +max
double getRand(double max, std::mt19937_64 & rng) {
    // Example:
    //  rng()/rng.max() is a number between 0 and 1
    //  multipled by two becomes a value between 0 and 2
    //  minus 1 makes it a random number between -1 and +1
    //  If max = 1.5 for example, the returned value is ( 1.5 * (value between -1 and +1) )
    return max * ( (static_cast<double>(rng()) / rng.max()) * 2.0 - 1);
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Creates a new rkParameters individual by combining properties of two parent Individuals using a mask to determine which
rkParameters<double> generateNewChild(const rkParameters<double> & p1, const rkParameters<double> & p2, const std::vector<int> & mask, const cudaConstants * cConstants, const double & annealing, std::mt19937_64 & rng, const double & generation) {
    // Set the new individual to hold traits from parent 1 
    rkParameters<double> newChild = p1;


    // Based on mask values, set parameters value is to be set to parent 2 or average of parent 1 and 2
    // Only calculate values pertaining to thrust if a thruster is being used
    if (cConstants->thruster_type != thruster<double>::NO_THRUST) {
        // Iterate through the gamma values
        for (int i = GAMMA_OFFSET; i < (GAMMA_OFFSET + GAMMA_ARRAY_SIZE); i++) {
            if (mask[i] == PARTNER2) {
                newChild.coeff.gamma[i - GAMMA_OFFSET] = p2.coeff.gamma[i - GAMMA_OFFSET]; // If the mask is 2, use the value from parent 2
            }
            else if (mask[i] == AVG) {
                // If the mask is 3, average the values from both parents
                newChild.coeff.gamma[i - GAMMA_OFFSET] = p2.coeff.gamma[i - GAMMA_OFFSET]/2.0 + p1.coeff.gamma[i - GAMMA_OFFSET]/2.0;
            }
        }
        // Iterate through tau values
        for (int i = TAU_OFFSET; i < (TAU_OFFSET + TAU_ARRAY_SIZE); i++) {
            if (mask[i] == PARTNER2) {
                newChild.coeff.tau[i - TAU_OFFSET] = p2.coeff.tau[i - TAU_OFFSET];
            }
            else if (mask[i] == AVG) {
                newChild.coeff.tau[i - TAU_OFFSET] = p2.coeff.tau[i - TAU_OFFSET]/2.0 + p1.coeff.tau[i - TAU_OFFSET]/2.0;
            }
        }
        // Iterate through coasting values
        for (int i = COAST_OFFSET; i < (COAST_OFFSET + COAST_ARRAY_SIZE); i++) {
            if (mask[i] == PARTNER2) {
                newChild.coeff.coast[i - COAST_OFFSET] = p2.coeff.coast[i - COAST_OFFSET];
            }
            else if (mask[i] == AVG) {
                newChild.coeff.coast[i - COAST_OFFSET] = p2.coeff.coast[i - COAST_OFFSET]/2.0 + p1.coeff.coast[i - COAST_OFFSET]/2.0;
            }
        }
    }
    
    // Go through other variables
    if (mask[TRIPTIME_OFFSET] == PARTNER2) { //tripTime
        newChild.tripTime = p2.tripTime;
    }
    else if (mask[TRIPTIME_OFFSET] == AVG) {
        newChild.tripTime = p2.tripTime/2.0 + p1.tripTime/2.0;
    }
    
    if (mask[ZETA_OFFSET] == PARTNER2) { //zeta
        newChild.zeta = p2.zeta;
    }
    else if (mask[ZETA_OFFSET] == AVG) {
        newChild.zeta = p1.zeta/2.0 + p2.zeta/2.0;
    }
    
    if (mask[BETA_OFFSET] == PARTNER2) { //beta
        newChild.beta = p2.beta;
    }
    else if (mask[BETA_OFFSET] == AVG) {
        newChild.beta = p2.beta/2.0 + p1.beta/2.0;
    }
    
    if (mask[ALPHA_OFFSET] == PARTNER2) { //alpha
        newChild.alpha = p2.alpha;
    }
    else if (mask[ALPHA_OFFSET] == AVG) {
        newChild.alpha = p1.alpha/2.0 + p2.alpha/2.0;
    }
    
    // Crossover complete, determine mutation
    newChild = mutate(newChild, rng, annealing, cConstants, generation);

    return newChild;    
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Utility function to generate a boolean mask that determines which parameter value is mutating and how many based on mutation_rate iteratively
void mutateMask(std::mt19937_64 & rng, bool * mutateMask, double mutation_rate) {
    for (int i = 0; i < OPTIM_VARS; i++) {
        //Reset mask
        mutateMask[i] = false;
    }
    
    // counter to make sure in the unlikely event that all genes are being mutated the code doesn't get stuck in infinite loop looking for a false spot in the array
    int geneCount = 0; 
    // Set a gene to mutate if a randomized values is less than mutation_rate, repeating everytime this is true
    while ((static_cast<double>(rng()) / rng.max()) < mutation_rate && geneCount < OPTIM_VARS) {
        bool geneSet = false; // boolean to flag if a gene was successfully selected to be set as true
        int index; // index value that will be assigned randomly to select a gene
        while (geneSet != true) {
            index = rng() % OPTIM_VARS;
            // If the randomly selected gene hasn't been already set to mutate, set it and flag the geneSet to true to break out of the loop
            if (mutateMask[index] == false) {
                mutateMask[index] = true;
                geneSet = true;
                geneCount++;
            }
        }

    }
}

// In a given Individual's parameters, generate a mutate mask using mutateMask() and then adjust parameters based on the mask, mutation of at least one gene is not guranteed
rkParameters<double> mutate(const rkParameters<double> & p1, std::mt19937_64 & rng, const double & annealing, const cudaConstants* cConstants, const double & generation) {    
    // initially set new individual to have all parameter values from parent 1
    rkParameters<double> newChild = p1;

    // Declare and set a mutation_mask for which gene is being mutated
    bool * mutation_mask = new bool[OPTIM_VARS];

    // Determine which parameters will be mutated with the mutate mask
    mutateMask(rng, mutation_mask, cConstants->mutation_rate);

    // Declare a record that is to describe what genes are being changed and by how much to record into mutateFile
    // double recordLog[OPTIM_VARS];

    // Iterate through the mutation_mask, mutating the corresponding gene if set to true
    for (int index = 0; index < OPTIM_VARS; index++) {
        if (mutation_mask[index] == true) {
            
            if ( (index >= GAMMA_OFFSET) && (index <= (GAMMA_OFFSET + GAMMA_ARRAY_SIZE-1)) ) { // Gamma value
                double randVar = getRand(cConstants->gamma_mutate_scale * annealing, rng);
                newChild.coeff.gamma[index-GAMMA_OFFSET] += randVar;
                // recordLog[index] = randVar;
            }
            else if ( (index >= TAU_OFFSET) && (index <= (TAU_OFFSET + TAU_ARRAY_SIZE-1))) { // Tau value 
                double randVar = getRand(cConstants->tau_mutate_scale * annealing, rng);
                newChild.coeff.tau[index-TAU_OFFSET] += randVar;
                // recordLog[index] = randVar;
            }
            else if (index >= COAST_OFFSET && index <= (COAST_OFFSET + COAST_ARRAY_SIZE-1)) { // Coast value
                double randVar = getRand(cConstants->coast_mutate_scale * annealing, rng);
                newChild.coeff.coast[index-COAST_OFFSET] += randVar;
                // recordLog[index] = randVar;
            }
            else if (index == TRIPTIME_OFFSET) { // Time final
                double randVar = getRand(cConstants->triptime_mutate_scale * annealing, rng);
                newChild.tripTime += randVar;
                // bound checking to make sure the tripTime is set within the valid range of trip times
                if (newChild.tripTime < cConstants->triptime_min + cConstants->timeRes) {
                    newChild.tripTime = cConstants->triptime_min + cConstants->timeRes;
                }
                else if (newChild.tripTime > cConstants->triptime_max - cConstants->timeRes) {
                    newChild.tripTime = cConstants->triptime_max - cConstants->timeRes;
                }

                // recordLog[index] = randVar;
            }
            else if (index == ZETA_OFFSET) { // Zeta
                double randVar = getRand(cConstants->zeta_mutate_scale * annealing, rng);
                newChild.zeta += randVar;
                // recordLog[index] = randVar;
            }
            else if (index == BETA_OFFSET) { // Beta
                double randVar = getRand(cConstants->beta_mutate_scale * annealing, rng);
                newChild.beta += randVar;
                // recordLog[index] = randVar;
    
                // A check to ensure beta remains in value range 0 to pi, doesn't update recordLog
                while (newChild.beta < 0 || newChild.beta > M_PI) {
                    if (newChild.beta < 0) {
                        newChild.beta += M_PI;
                    }
                    else if (newChild.beta > M_PI) {
                        newChild.beta -= M_PI;              
                    }
                }
            }
            else if (index == ALPHA_OFFSET) { // Alpha
                double randVar = getRand(cConstants->alpha_mutate_scale * annealing, rng);
                newChild.alpha += randVar;                
                // recordLog[index] = randVar;
            }
        }
        // else { // Record if the gene is not being mutated
        //     recordLog[index] = 0;
        // }
    }

    // If in record mode, append the recordLog into the .csv file
    // if (cConstants->record_mode == true) {
    //     int genesMutated = 0;
    //     for (int i = 0; i < OPTIM_VARS; i++) {
    //         if (mutation_mask[i] == true) {
    //             genesMutated++;
    //         }
    //     }
    //     recordMutateFile(cConstants, generation, annealing, genesMutated, recordLog);
    // }

    //Delete the mask and set of parameters to free memory
    delete [] mutation_mask;
    return newChild;
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Method that creates a pair of new Individuals from a pair of other individuals and a mask
void generateChildrenPair (Adult parent1, Adult parent2, Child * newChildren, std::vector<int> & mask, const double & annealing, std::mt19937_64 & rng, int & numNewChildren, const int & generation, const cudaConstants* cConstants) {
    //Generate a new individual based on the two parents
    newChildren[numNewChildren] = Child(generateNewChild(parent1.startParams, parent2.startParams, mask,  cConstants, annealing, rng, generation), cConstants);

    //Add one to numNewChildren to signify that a new child has been added to the newChildren vector
    //Will also ensure that no children are overwritten
    numNewChildren++;

    //Flip the mask to generate a mirrored individual
    flipMask(mask); 

    //Generate a mirrored child
    newChildren[numNewChildren] = Child(generateNewChild(parent1.startParams, parent2.startParams, mask,  cConstants, annealing, rng, generation), cConstants); 

    //Add one to numNewChildren since a new child was generated
    numNewChildren++; 
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Creates the next pool to be used in the optimize function in opimization.cu
void newGeneration (std::vector<Adult> & oldAdults, std::vector<Adult> & newAdults, const double & annealing, const int & generation, std::mt19937_64 & rng, const cudaConstants* cConstants) {
    //Import number of new individuals the GPU needs to fill its threads
    //int newChildren = cConstants -> num_individuals;

    //Create a newChildren function to fill in with children generated from oldAdults
    //Child * newChildren[newChildren];

    //Import number of new individuals the GPU needs to fill its threads
    //Create a newChildren function to fill in with children generated from oldAdults
    Child* newChildren = new Child[cConstants->num_individuals]; 

    //Create a mask to determine which of the child's parameters will be inherited from which parents
    std::vector<int> mask;

    //Shuffle the oldAdults vector to randomize the parents for each of the children
    std::shuffle(std::begin(oldAdults), std::end(oldAdults), rng); //possibly need to use iterators for the beginning and ending and a timeseed (std::default_random_engine(seed)) -> https://www.cplusplus.com/reference/algorithm/shuffle/?kw=shuffle

    //Create a int that will store the number of pairs needed to generate enough children for the number of threads
    //Assigned as num_individuals / 8 as each pair will generate 8 individuals
    int numPairs = cConstants->num_individuals / 8;

    //Int that tracks the number of new individuals that have been generated
    //Used primarily to make sure that no children are overwritten in newChildren; indexing in generatechildrenPair, should equal N when finished
    //Starts as 0 beacuse no new children have been generated
    int numNewChildren = 0;

    //Generate a pair of children based on the random cross over mask method for each pair of parents
    for (int i = 0; i < numPairs; i += 2) {
        //Generate the base mask
        crossOver_wholeRandom(mask, rng);

        //Generate a pair of children based on the mask
        generateChildrenPair(oldAdults[i], oldAdults[i+1], newChildren, mask, annealing, rng, numNewChildren, generation, cConstants);
    }
    
    //Generate a pair of children for each pair of parents based on the average mask method
    for (int i = 0; i < numPairs; i += 2) {
        //Generate the base mask
        crossOver_average(mask);

        //Generate a pair of children based on the mask
        generateChildrenPair(oldAdults[i], oldAdults[i+1], newChildren, mask, annealing, rng, numNewChildren, generation, cConstants);
    }

    //Generate a pair of children for each pair of parents based on the bundled random mask method
    for (int i = 0; i < numPairs; i += 2) {
        //Generate the base mask
        crossOver_bundleVars(mask, rng);

        //Generate a pair of children based on the mask
        generateChildrenPair(oldAdults[i], oldAdults[i+1], newChildren, mask, annealing, rng, numNewChildren, generation, cConstants);
    }
    
    //Generate a pair of children for each pair of parents based on the bundled random mask method for a second time
    for (int i = 0; i < numPairs; i += 2) {
        //Generate the base mask
        crossOver_bundleVars(mask, rng);

        //Generate a pair of children based on the mask
        generateChildrenPair(oldAdults[i], oldAdults[i+1], newChildren, mask, annealing, rng, numNewChildren, generation, cConstants);
    }

    double timeInitial = 0;
    double calcPerS = 0;

    // Each child represents an individual set of starting parameters 
        // GPU based runge kutta process determines final position and velocity based on parameters
    //Will fill in the final variables (final position & speed, posDiff, speedDiff) for each child
    //TODO: get rid of magic numbers, first 0 is timeInitial, second 0 is calcPerS, the are both set as 0 before being passed in to callRK within optimize.cu
        //Perhaps we should just import cCOnstants and newChildren into callRk, since most of the arguments come from cConstants regardless
    callRK(cConstants->num_individuals, cConstants->thread_block_size, newChildren, timeInitial, (cConstants->orbitalPeriod / cConstants->GuessMaxPossibleSteps), cConstants->rk_tol, calcPerS, cConstants); 

    //Now that the children have been simulated, convert the children into adults
    //This will also put the converted children into the newAdults vector
    convertToAdults(newAdults, newChildren, cConstants); 

    //Free the pointer's memory
    delete[] newChildren; 
}

//Function that will convert the generated children into adults
//Will also transfer the created adults into the newAdult vector
void convertToAdults(std::vector<Adult> & newAdults, Child* newChildren, const cudaConstants* cConstants) {
    //Iterate through the newChildren vector to add them to the newAdult vector
    //iterates until num_individuals as that is the size of newChildren
    for (int i = 0; i < cConstants->num_individuals; i++)
    {
        //Fill in the newAdults vector with a new adult generated with a completed child
        newAdults.push_back(Adult(newChildren[i]));
    }
}

//Will create the first generation of adults from random parameters so that future generations have a pre-existing base of adults
void firstGeneration(Child* initialChildren, std::vector<Adult>& oldAdults, const cudaConstants* cConstants){
    const double timeIntial = 0;
    double calcPerS  = 0;
     // Each child represents an individual set of starting parameters 
        // GPU based runge kutta process determines final position and velocity based on parameters
    //Will fill in the final variables (final position & speed, posDiff, speedDiff) for each child
    //TODO: get rid of magic numbers, first 0 is timeInitial, second 0 is calcPerS, the are both set as 0 before being passed in to callRK within optimize.cu
        //Perhaps we should just import cCOnstants and newChildren into callRk, since most of the arguments come from cConstants regardless
    callRK(cConstants->num_individuals, cConstants->thread_block_size, initialChildren, timeIntial, (cConstants->orbitalPeriod / cConstants->GuessMaxPossibleSteps), cConstants->rk_tol, calcPerS, cConstants); 

    //Now that the children have been simulated, convert the children into adults
    //This will also put the converted children into the newAdults vector
    convertToAdults(oldAdults, initialChildren, cConstants); 
}