#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>

// Global enumeration for the different mask values instead of 1,2,3 for better readibility and clarity of value meaning
enum maskValue {
    PARTNER1 = 1,
    PARTNER2,
    AVG,
};

// Determing selection of survivors that will carry properties into the new individuals of the newGeneration
void selectSurvivors(Individual * pool, int poolSize, int selectionSize, Individual* survivors, const double & ratio, const int & missionType) {
    // Sort the pool by position difference
    // and assign the first part of the survivor array for best posDiff individuals
    // portion size based on the ratio percentage
    std::sort(pool, pool+poolSize, LowerPosDiff);
    //Select survivors (starting at 0)
    for (int i = 0; i < selectionSize*ratio; i++) {
        survivors[i] = pool[i];
    }

    // Sort the pool by speed difference. If mission type is soft, use LowerSpeedDiff because we want velocity = 0
    if(missionType == Impact){
        std::sort(pool, pool+poolSize, HigherSpeedDiff);
    }
    else if(missionType == Rendezvous){
        std::sort(pool, pool+poolSize, LowerSpeedDiff);
    }

    //Used to make sure pool[] starts at 0 
    int j = selectionSize*ratio; 
    //starting where first loop ended
    for (int i = selectionSize*ratio; i < selectionSize; i++) {
        survivors[(i)] = pool[i - j];
    }
    return;
}

// Creates a random bifurcation mask, currently not in use
void crossOver_randHalf(int * mask, std::mt19937_64 & rng) {
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
void crossOver_oneParent(int * mask) {
    for (int i = 0; i < OPTIM_VARS; i++) {
        mask[i] = PARTNER1;
    }
    return;
}

// Creates a random mask
void crossOver_wholeRandom(int * mask, std::mt19937_64 & rng) {
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
void crossOver_bundleVars(int * mask, std::mt19937_64 & rng) {
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
void crossOver_average(int * mask) {
    for (int i = 0; i < OPTIM_VARS; i++) {
        mask[i] = AVG;
    }
    return;
}

// Utility to flip the polarity of a mask
void flipMask(int * mask) {
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

// Creates a new rkParameters individual by combining properties of two parent Individuals using a mask to determine which
rkParameters<double> generateNewIndividual(const rkParameters<double> & p1, const rkParameters<double> & p2, const int * mask, const cudaConstants * cConstants, double annealing, std::mt19937_64 & rng, double generation, int gene) {
    // Set the new individual to hold traits from parent 1 
    rkParameters<double> newInd = p1;


    // Based on mask values, set parameters value is to be set to parent 2 or average of parent 1 and 2
    // Only calculate values pertaining to thrust if a thruster is being used
    if (cConstants->thruster_type != thruster<double>::NO_THRUST) {
        // Iterate through the gamma values
        for (int i = GAMMA_OFFSET; i < (GAMMA_OFFSET + GAMMA_ARRAY_SIZE); i++) {
            if (mask[i] == PARTNER2) {
                newInd.coeff.gamma[i - GAMMA_OFFSET] = p2.coeff.gamma[i - GAMMA_OFFSET]; // If the mask is 2, use the value from parent 2
            }
            else if (mask[i] == AVG) {
                // If the mask is 3, average the values from both parents
                newInd.coeff.gamma[i - GAMMA_OFFSET] = p2.coeff.gamma[i - GAMMA_OFFSET]/2.0 + p1.coeff.gamma[i - GAMMA_OFFSET]/2.0;
            }
        }
        // Iterate through tau values
        for (int i = TAU_OFFSET; i < (TAU_OFFSET + TAU_ARRAY_SIZE); i++) {
            if (mask[i] == PARTNER2) {
                newInd.coeff.tau[i - TAU_OFFSET] = p2.coeff.tau[i - TAU_OFFSET];
            }
            else if (mask[i] == AVG) {
                newInd.coeff.tau[i - TAU_OFFSET] = p2.coeff.tau[i - TAU_OFFSET]/2.0 + p1.coeff.tau[i - TAU_OFFSET]/2.0;
            }
        }
        // Iterate through coasting values
        for (int i = COAST_OFFSET; i < (COAST_OFFSET + COAST_ARRAY_SIZE); i++) {
            if (mask[i] == PARTNER2) {
                newInd.coeff.coast[i - COAST_OFFSET] = p2.coeff.coast[i - COAST_OFFSET];
            }
            else if (mask[i] == AVG) {
                newInd.coeff.coast[i - COAST_OFFSET] = p2.coeff.coast[i - COAST_OFFSET]/2.0 + p1.coeff.coast[i - COAST_OFFSET]/2.0;
            }
        }
    }
    
    // Go through other variables
    if (mask[TRIPTIME_OFFSET] == PARTNER2) { //tripTime
        newInd.tripTime = p2.tripTime;
    }
    else if (mask[TRIPTIME_OFFSET] == AVG) {
        newInd.tripTime = p2.tripTime/2.0 + p1.tripTime/2.0;
    }
    
    if (mask[ZETA_OFFSET] == PARTNER2) { //zeta
        newInd.zeta = p2.zeta;
    }
    else if (mask[ZETA_OFFSET] == AVG) {
        newInd.zeta = p1.zeta/2.0 + p2.zeta/2.0;
    }
    
    if (mask[BETA_OFFSET] == PARTNER2) { //beta
        newInd.beta = p2.beta;
    }
    else if (mask[BETA_OFFSET] == AVG) {
        newInd.beta = p2.beta/2.0 + p1.beta/2.0;
    }
    
    if (mask[ALPHA_OFFSET] == PARTNER2) { //alpha
        newInd.alpha = p2.alpha;
    }
    else if (mask[ALPHA_OFFSET] == AVG) {
        newInd.alpha = p1.alpha/2.0 + p2.alpha/2.0;
    }
    
    // Crossover complete, determine mutation
    newInd = mutate(newInd, rng, annealing, cConstants, generation, gene);

    return newInd;    
}

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

// Utility function to generate a boolean mask that determines which parameter value is mutating and how many based on mutation_rate iteratively
void mutateGeneMask(std::mt19937_64 & rng, bool * mutateMask, double mutation_rate, int gene) {
    for (int i = 0; i < OPTIM_VARS; i++) {
        //Reset mask
        mutateMask[i] = false;
    }
    if (gene == GAMMA_OFFSET){
        for(int i = GAMMA_OFFSET; i < GAMMA_OFFSET + GAMMA_ARRAY_SIZE; i++){
            mutateMask[i] = true;
        }
    }
    else if (gene == TAU_OFFSET){
        for(int i = TAU_OFFSET; i < TAU_OFFSET + TAU_ARRAY_SIZE; i++){
            mutateMask[i] = true;
        }
    }
    else if (gene == COAST_OFFSET){
        for(int i = COAST_OFFSET; i < COAST_OFFSET + COAST_ARRAY_SIZE; i++){
            mutateMask[i] = true;
        }
    }
    else {
        //Set requested gene to be mutated 
        mutateMask[gene] = true;
    }

}

// In a given Individual's parameters, generate a mutate mask using mutateMask() and then adjust parameters based on the mask, mutation of at least one gene is not guranteed
rkParameters<double> mutate(const rkParameters<double> & p1, std::mt19937_64 & rng, double annealing, const cudaConstants* cConstants, double generation, int gene) {    
    // initially set new individual to have all parameter values from parent 1
    rkParameters<double> newInd = p1;

    // Declare and set a mutation_mask for which gene is being mutated
    bool * mutation_mask = new bool[OPTIM_VARS];
    if(gene == -1){
        mutateMask(rng, mutation_mask, cConstants->mutation_rate);
    }
    else{
        mutateGeneMask(rng, mutation_mask, cConstants->mutation_rate, gene);
    }


    // Declare a record that is to describe what genes are being changed and by how much to record into mutateFile
    // double recordLog[OPTIM_VARS];

    // Iterate through the mutation_mask, mutating the corresponding gene if set to true
    for (int index = 0; index < OPTIM_VARS; index++) {
        if (mutation_mask[index] == true) {
            
            if ( (index >= GAMMA_OFFSET) && (index <= (GAMMA_OFFSET + GAMMA_ARRAY_SIZE-1)) ) { // Gamma value
                double randVar = getRand(cConstants->gamma_mutate_scale * annealing, rng);
                newInd.coeff.gamma[index-GAMMA_OFFSET] += randVar;
                // recordLog[index] = randVar;
            }
            else if ( (index >= TAU_OFFSET) && (index <= (TAU_OFFSET + TAU_ARRAY_SIZE-1))) { // Tau value 
                double randVar = getRand(cConstants->tau_mutate_scale * annealing, rng);
                newInd.coeff.tau[index-TAU_OFFSET] += randVar;
                // recordLog[index] = randVar;
            }
            else if (index >= COAST_OFFSET && index <= (COAST_OFFSET + COAST_ARRAY_SIZE-1)) { // Coast value
                double randVar = getRand(cConstants->coast_mutate_scale * annealing, rng);
                newInd.coeff.coast[index-COAST_OFFSET] += randVar;
                // recordLog[index] = randVar;
            }
            else if (index == TRIPTIME_OFFSET) { // Time final
                double randVar = getRand(cConstants->triptime_mutate_scale * annealing, rng);
                newInd.tripTime += randVar;
                // bound checking to make sure the tripTime is set within the valid range of trip times
                if (newInd.tripTime < cConstants->triptime_min + cConstants->timeRes) {
                    newInd.tripTime = cConstants->triptime_min + cConstants->timeRes;
                }
                else if (newInd.tripTime > cConstants->triptime_max - cConstants->timeRes) {
                    newInd.tripTime = cConstants->triptime_max - cConstants->timeRes;
                }

                // recordLog[index] = randVar;
            }
            else if (index == ZETA_OFFSET) { // Zeta
                double randVar = getRand(cConstants->zeta_mutate_scale * annealing, rng);
                newInd.zeta += randVar;
                // recordLog[index] = randVar;
            }
            else if (index == BETA_OFFSET) { // Beta
                double randVar = getRand(cConstants->beta_mutate_scale * annealing, rng);
                newInd.beta += randVar;
                // recordLog[index] = randVar;
    
                // A check to ensure beta remains in value range 0 to pi, doesn't update recordLog
                if (newInd.beta < 0) {
                    newInd.beta = 0;
                }
                else if (newInd.beta > M_PI) {
                    newInd.beta = M_PI;              
                }
            }
            else if (index == ALPHA_OFFSET) { // Alpha
                double randVar = getRand(cConstants->alpha_mutate_scale * annealing, rng);
                newInd.alpha += randVar;                
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

    delete [] mutation_mask;
    return newInd;
}


// Method that creates a pair of new Individuals from a pair of other individuals and a mask
void generateChildrenPair(Individual *pool, Individual *survivors, int * mask, int& newIndCount, int parentsIndex, double annealing, int poolSize, std::mt19937_64 & rng, const cudaConstants* cConstants, double generation, int gene) { 
    // Determine where the parents and the new individual being created are located in the pool
    int parent1Index = parentsIndex;
    int parent2Index = parentsIndex + 1;

    // The new indiviudal is located at the end of the pool
    // up the number of new individuals already created
    int newIndividualIndex = poolSize - 1 - newIndCount;    
    // Generate new offspring with mask
    pool[newIndividualIndex] = Individual(generateNewIndividual(survivors[parent1Index].startParams, survivors[parent2Index].startParams, mask, cConstants, annealing, rng, generation, gene), cConstants);
    newIndCount++;

    // Get the opposite offspring from the mask by flipping the mask
    newIndividualIndex--; // Decrement newIndividualIndex value to access where the next individual must be as newIndCount has increased
    flipMask(mask);
    pool[newIndividualIndex] = Individual(generateNewIndividual(survivors[parent1Index].startParams, survivors[parent2Index].startParams, mask, cConstants, annealing, rng, generation, gene), cConstants);
    newIndCount++;

    return;
}

// Creates the next pool to be used in the optimize function in opimization.cu
int newGeneration(Individual *survivors, Individual *pool, int survivorSize, int poolSize, double annealing, const cudaConstants* cConstants, std::mt19937_64 & rng, double generation) {
    //Crossover mask allocation, used for performance
    int * mask = new int[OPTIM_VARS];
    // Number of new individuals created so far (initially none)
    // used in navigating through the pool when creating new individuals and returned at end of function
    int newIndCount = 0;

    //-------------------------------------NEW METHOD--------------------------------------------------------------------------------------------------
    // If a generation gets a certain amount of clones, use this crossover method instead of the normal one to replace the clones.
    // For this method to work, it is important that the pool is sorted by cost. That should be happening in optimization because it
    // is also important for the normal crossovers.

    //number of clones that can exist before needing to be replaced.
    int cloneThresh = 1439;

    //count of clones identical to the best cost clone.
    int cloneCount = 0;
     
    // Loop through all individuals. If they are a clone, increase cloneCount and toggle isClone bool for that individual.
    for(int i = 1; i < poolSize; i++){
        if (same(pool[0], pool[i], 0.01)){
            cloneCount++;
            pool[i].isClone = true;
        }
    }

    std::cout << cloneCount;
    // If cloneCount reaches a certain threshold, replace the clones with mutated versions of themselves.
    // if (cloneCount > cloneThresh){
    //     //Set mask to do no crossover and copy the same individual
    //     crossOver_oneParent(mask);
    //     //loop through all the individuals and find the ones that are clones
    //     for (int i=1; i < poolSize; i++){
    //         if(pool[i].isClone == true){
    //             //replace individual with mutated version of itself.
    //             int index = rng() % 7;
    //             switch (index){
    //                 case '0':
    //                     pool[i] = Individual(generateNewIndividual(pool[i].startParams, pool[i].startParams, mask, cConstants, annealing, rng, generation, GAMMA_OFFSET), cConstants);
    //                     break;
    //                 case '1':
    //                     pool[i] = Individual(generateNewIndividual(pool[i].startParams, pool[i].startParams, mask, cConstants, annealing, rng, generation, TAU_OFFSET), cConstants);
    //                     break;
    //                 case '2':
    //                     pool[i] = Individual(generateNewIndividual(pool[i].startParams, pool[i].startParams, mask, cConstants, annealing, rng, generation, ALPHA_OFFSET), cConstants);
    //                     break;
    //                 case '3':
    //                     pool[i] = Individual(generateNewIndividual(pool[i].startParams, pool[i].startParams, mask, cConstants, annealing, rng, generation, BETA_OFFSET), cConstants);
    //                     break;
    //                 case '4':
    //                     pool[i] = Individual(generateNewIndividual(pool[i].startParams, pool[i].startParams, mask, cConstants, annealing, rng, generation, ZETA_OFFSET), cConstants);
    //                     break;
    //                 case '5':
    //                     pool[i] = Individual(generateNewIndividual(pool[i].startParams, pool[i].startParams, mask, cConstants, annealing, rng, generation, TRIPTIME_OFFSET), cConstants);
    //                     break;
    //                 case '6':
    //                     pool[i] = Individual(generateNewIndividual(pool[i].startParams, pool[i].startParams, mask, cConstants, annealing, rng, generation, COAST_OFFSET), cConstants);
    //                     break;
    //             }
    //             //increment new individual count
    //             newIndCount++; 
    //         }
    //     }
    //     //Sort based on cost. This will put the new individuals at the bottom of the pack so that they get replaced. (The constructor's cost = 10)
    //     std::sort(pool, pool + poolSize);
    // }
    //else {//Normal crossover methods when clones have not exeeded cloneThresh
    
        // Value for how many pairs to use and produce in each loop
        // (as one iteration through a loop produces a new pair)
        int numPairs = survivorSize / 2;

        // Shuffle the survivors to ensure diverse crossover
        std::shuffle(survivors, survivors+survivorSize, rng);

        // Generate two offspring through each crossover method
        // total is 4 * survivorSize offspring in pool

        // Every loop needs a unique mask each iteration

        // Loop for wholeRandom mask
        for (int i = 0; i < numPairs; i++) {
            crossOver_wholeRandom(mask, rng);
            generateChildrenPair(pool, survivors, mask, newIndCount, 2*i, annealing, poolSize, rng, cConstants, generation, -1);
        }
        // Loop for averaging mask
        for (int i = 0; i < numPairs; i++) {
            crossOver_average(mask);
            generateChildrenPair(pool, survivors, mask, newIndCount, 2*i, annealing, poolSize, rng, cConstants, generation, -1);
        }
        // 2 loops for bundleVars mask,
        // two seperate loops resulting from carry over of past code
        // also will allow easier changes in masks used (right now just using two bundleVars instead of two different ones)
        for (int i = 0; i < numPairs; i++) {
            crossOver_bundleVars(mask, rng);
            generateChildrenPair(pool, survivors, mask, newIndCount, 2*i, annealing, poolSize, rng, cConstants, generation, -1);
        }
        for (int i = 0; i < numPairs; i++) {
            crossOver_bundleVars(mask, rng);
            generateChildrenPair(pool, survivors, mask, newIndCount, 2*i, annealing, poolSize, rng, cConstants, generation, -1);
        }
    //}

    delete [] mask;
    return newIndCount;
}
