#include <algorithm>
#include <iostream>

// Global enumeration for the different mask values instead of 1,2,3 for better readibility and clarity of value meaning
enum maskValue {
    PARTNER1 = 1,
    PARTNER2,
    AVG,
    AVG_RATIO,
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

// Sets the entire mask to be AVG for length OPTIM_VARS
void crossOver_averageRatio(std::vector<int> & mask) {
    for (int i = 0; i < OPTIM_VARS; i++) {
        mask[i] = AVG_RATIO;
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
        else if (mask[i] == AVG){
            mask[i] = AVG_RATIO;
        }
        else if (mask[i] == AVG_RATIO){
            mask[i] = AVG;
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
rkParameters<double> generateNewChild(const rkParameters<double> & p1, const rkParameters<double> & p2, const std::vector<int> & mask, const cudaConstants * cConstants, const double & annealing, std::mt19937_64 & rng, const int & generation) {

    // Set the new individual to hold traits from parent 1 
    rkParameters<double> childParameters = p1;

    double ratio;
    ratio = (static_cast<double>(rng()) / rng.max());//random value between 0 and 1
    // Based on mask values, set parameters value is to be set to parent 2 or average of parent 1 and 2
    // Only calculate values pertaining to thrust if a thruster is being used
    if (cConstants->thruster_type != thruster<double>::NO_THRUST) {
        // Iterate through the gamma values
        for (int i = GAMMA_OFFSET; i < (GAMMA_OFFSET + GAMMA_ARRAY_SIZE); i++) {
            if (mask[i] == PARTNER2) {
                childParameters.coeff.gamma[i - GAMMA_OFFSET] = p2.coeff.gamma[i - GAMMA_OFFSET]; // If the mask is 2, use the value from parent 2
            }
            else if (mask[i] == AVG) {
                // If the mask is 3, average the values from both parents
                childParameters.coeff.gamma[i - GAMMA_OFFSET] = p2.coeff.gamma[i - GAMMA_OFFSET]/2.0 + p1.coeff.gamma[i - GAMMA_OFFSET]/2.0;
            }
            else if (mask[i] == AVG_RATIO){
                // If the mask is 4, average the values from both parents based on random ratio
                childParameters.coeff.gamma[i - GAMMA_OFFSET] = ratio*p2.coeff.gamma[i - GAMMA_OFFSET] + (1 - ratio)*p1.coeff.gamma[i - GAMMA_OFFSET];
            }
        }
        // Iterate through tau values
        for (int i = TAU_OFFSET; i < (TAU_OFFSET + TAU_ARRAY_SIZE); i++) {
            if (mask[i] == PARTNER2) {
                childParameters.coeff.tau[i - TAU_OFFSET] = p2.coeff.tau[i - TAU_OFFSET];
            }
            else if (mask[i] == AVG) {
                childParameters.coeff.tau[i - TAU_OFFSET] = p2.coeff.tau[i - TAU_OFFSET]/2.0 + p1.coeff.tau[i - TAU_OFFSET]/2.0;
            }
            else if (mask[i] == AVG_RATIO) {
                childParameters.coeff.tau[i - TAU_OFFSET] = ratio*p2.coeff.tau[i - TAU_OFFSET] + (1 - ratio)*p1.coeff.tau[i - TAU_OFFSET];
            }
        }
        // Iterate through coasting values
        for (int i = COAST_OFFSET; i < (COAST_OFFSET + COAST_ARRAY_SIZE); i++) {
            if (mask[i] == PARTNER2) {
                childParameters.coeff.coast[i - COAST_OFFSET] = p2.coeff.coast[i - COAST_OFFSET];
            }
            else if (mask[i] == AVG) {
                childParameters.coeff.coast[i - COAST_OFFSET] = p2.coeff.coast[i - COAST_OFFSET]/2.0 + p1.coeff.coast[i - COAST_OFFSET]/2.0;
            }
            else if (mask[i] == AVG_RATIO) {
                childParameters.coeff.coast[i - COAST_OFFSET] = ratio*p2.coeff.coast[i - COAST_OFFSET] + (1 - ratio)*p1.coeff.coast[i - COAST_OFFSET];
            }
        }
    }
    
    // Go through other variables
    if (mask[TRIPTIME_OFFSET] == PARTNER2) { //tripTime
        childParameters.tripTime = p2.tripTime;
    }
    else if (mask[TRIPTIME_OFFSET] == AVG) {
        childParameters.tripTime = p2.tripTime/2.0 + p1.tripTime/2.0;
    }
    else if (mask[TRIPTIME_OFFSET] == AVG_RATIO) {
        childParameters.tripTime = ratio*p2.tripTime + (1 - ratio)*p1.tripTime;
    }
    
    if (mask[ZETA_OFFSET] == PARTNER2) { //zeta
        childParameters.zeta = p2.zeta;
    }
    else if (mask[ZETA_OFFSET] == AVG) {
        childParameters.zeta = p2.zeta/2.0 + p1.zeta/2.0;
    }
    else if (mask[ZETA_OFFSET] == AVG_RATIO) {
        childParameters.zeta = ratio*p2.zeta + (1 - ratio)*p1.zeta;
    }
    
    if (mask[BETA_OFFSET] == PARTNER2) { //beta
        childParameters.beta = p2.beta;
    }
    else if (mask[BETA_OFFSET] == AVG) {
        childParameters.beta = p2.beta/2.0 + p1.beta/2.0;
    }
    else if (mask[BETA_OFFSET] == AVG_RATIO) {
        childParameters.beta = ratio*p2.beta + (1 - ratio)*p1.beta;
    }
    
    if (mask[ALPHA_OFFSET] == PARTNER2) { //alpha
        childParameters.alpha = p2.alpha;
    }
    else if (mask[ALPHA_OFFSET] == AVG) {
        childParameters.alpha = p2.alpha/2.0 + p1.alpha/2.0;
    }
    else if (mask[ALPHA_OFFSET] == AVG_RATIO) {
        childParameters.alpha = ratio*p2.alpha + (1 - ratio)*p1.alpha;
    }

    // Crossover complete, determine mutation
    childParameters = mutate(childParameters, rng, annealing, cConstants, generation, cConstants->mutation_amplitude, cConstants->default_mutation_chance);

    return childParameters;
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Utility function to generate a boolean mask that determines which parameter value is mutating and how many based on mutation_rate iteratively
void mutateMask(std::mt19937_64 & rng, bool * mutateMask, double mutation_rate) {
    for (int i = 0; i < OPTIM_VARS; i++) {
        //Reset mask
        mutateMask[i] = false;
    }
    
    // counter to make sure in the unlikely event that all genes are being mutated the code doesn't
    // get stuck in infinite loop looking for a false spot in the array
    int geneCount = 0; //How many gene's have we changed
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

//TODO: take out duplicate
//      use two different parameters for both mutation chance (we probably want higher for duplicate mutiation)
//                                        and mutationScale (the amount we mutate)

// In a given Individual's parameters, generate a mutate mask using mutateMask() and then adjust parameters based on the mask, mutation of at least one gene is not guranteed
rkParameters<double> mutate(const rkParameters<double> & p1, std::mt19937_64 & rng, const double & annealing, const cudaConstants* cConstants, const int & generation, const double & mutationScale, const double & mutation_chance) {    
    // initially set new individual to have all parameter values from parent 1
    rkParameters<double> childParameters = p1; 

    // Declare and set a mutation_mask for which gene is being mutated
    bool * mutation_mask = new bool[OPTIM_VARS];

    // Determine which parameters will be mutated with the mutate mask
    mutateMask(rng, mutation_mask, mutation_chance);
    

    // Declare a record that is to describe what genes are being changed and by how much to record into mutateFile
    // double recordLog[OPTIM_VARS];

    // Iterate through the mutation_mask, mutating the corresponding gene if set to true
    for (int index = 0; index < OPTIM_VARS; index++) {
        if (mutation_mask[index] == true) {
            
            if ( (index >= GAMMA_OFFSET) && (index <= (GAMMA_OFFSET + GAMMA_ARRAY_SIZE-1)) ) { // Gamma value
                double randVar = getRand(cConstants->gamma_mutate_scale * annealing * mutationScale, rng);
                childParameters.coeff.gamma[index-GAMMA_OFFSET] += randVar;
                // recordLog[index] = randVar;
            }
            else if ( (index >= TAU_OFFSET) && (index <= (TAU_OFFSET + TAU_ARRAY_SIZE-1))) { // Tau value 
                double randVar = getRand(cConstants->tau_mutate_scale * annealing * mutationScale, rng);
                childParameters.coeff.tau[index-TAU_OFFSET] += randVar;
                // recordLog[index] = randVar;
            }
            else if (index >= COAST_OFFSET && index <= (COAST_OFFSET + COAST_ARRAY_SIZE-1)) { // Coast value
                double randVar = getRand(cConstants->coast_mutate_scale * annealing * mutationScale, rng);
                childParameters.coeff.coast[index-COAST_OFFSET] += randVar;
                // recordLog[index] = randVar;
            }
            else if (index == TRIPTIME_OFFSET) { // Time final
                double randVar = getRand(cConstants->triptime_mutate_scale * annealing * mutationScale, rng);
                childParameters.tripTime += randVar;
                // bound checking to make sure the tripTime is set within the valid range of trip times
                if (childParameters.tripTime < cConstants->triptime_min + cConstants->timeRes) {
                    childParameters.tripTime = cConstants->triptime_min + cConstants->timeRes;
                }
                else if (childParameters.tripTime > cConstants->triptime_max - cConstants->timeRes) {
                    childParameters.tripTime = cConstants->triptime_max - cConstants->timeRes;
                }

                // recordLog[index] = randVar;
            }
            else if (index == ZETA_OFFSET) { // Zeta
                double randVar = getRand(cConstants->zeta_mutate_scale * annealing * mutationScale, rng);
                childParameters.zeta += randVar;
                // recordLog[index] = randVar;
            }
            else if (index == BETA_OFFSET) { // Beta
                double randVar = getRand(cConstants->beta_mutate_scale * annealing * mutationScale, rng);
                childParameters.beta += randVar;
                // recordLog[index] = randVar;
    
                // A check to ensure beta remains in value range 0 to pi, doesn't update recordLog
                while (childParameters.beta < 0 || childParameters.beta > M_PI) {
                    if (childParameters.beta < 0) {
                        childParameters.beta += M_PI;
                    }
                    else if (childParameters.beta > M_PI) {
                        childParameters.beta -= M_PI;              
                    }
                }
            }
            else if (index == ALPHA_OFFSET) { // Alpha
                double randVar = getRand(cConstants->alpha_mutate_scale * annealing * mutationScale, rng);
                childParameters.alpha += randVar;                
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
    return childParameters;
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Method that creates a pair of new Individuals from a pair of other individuals and a mask
void generateChildrenPair (const Adult & parent1, const Adult & parent2, Child * newChildren, const int & childrenToGenerate, std::vector<int> & mask, const double & annealing, std::mt19937_64 & rng, int & numNewChildren, const int & generation, const cudaConstants* cConstants) {
    
    //checks to make sure that more children are needed in the newChildren array before generating new children
    //if there are already enough children, it just returns and exits this function
    if(numNewChildren >= childrenToGenerate){ 
        return;
    }


    //this is not a duplicate

    //TODO: There is inherently an issue here that you may over-write newChildren if you don't know the size
    //      It is possible you can do this externally, but here in this function, there is no check on the size of numNewChildren
    //Generate a new individual based on the two parents
    newChildren[numNewChildren] = Child(generateNewChild(parent1.startParams, parent2.startParams, mask, cConstants, annealing, rng, generation), cConstants, generation);

    //Add one to numNewChildren to signify that a new child has been added to the newChildren vector
    //Will also ensure that no children are overwritten
    numNewChildren++;

    //checks to make sure that more children are needed in the newChildren array before generating new children
    //if there are already enough children, it just returns and exits this function
    if(numNewChildren >= childrenToGenerate){
        return;
    }

    //TODO: For efficiency, we might want to check if the mask is an averaging mask (if it is, we probably don't need to make the 'flip')
    //      On the other hand, even numbers are nice, so we may want to take the hit anyway just to keep things even
    //Dealt with this by "flipping" an average mask to be a weighted average mask
    //Flip the mask to generate a mirrored individual
    flipMask(mask); 

    //Generate a mirrored child
    newChildren[numNewChildren] = Child(generateNewChild(parent1.startParams, parent2.startParams, mask, cConstants, annealing, rng, generation), cConstants, generation); 

    //Add one to numNewChildren since a new child was generated
    numNewChildren++;
}

// Generate children from non-duplicate parents using crossover methods
void generateChildrenFromCrossover(std::vector<Adult> &parents, Child *newChildren, const int &childrenToGenerate, std::mt19937_64 &rng, const double &currentAnneal, const int &generation, const cudaConstants *cConstants)
{
    //std::cout << "\n_-_-_-_-_-_-_-_-_-Start of generate children from crossover_-_-_-_-_-_-_-_-_-\n";
    //std::cout << "\n_-_-_-_-_-_-_-_-_-Size of parents: " << parents.size() << ", number of children to generate: " << childrenToGenerate << "_-_-_-_-_-_-_-_-_-\n";
    // Create a mask to determine which of the child's parameters will be inherited from which parents
    std::vector<int> mask;

    // Make a basic mask with it set to avg at first
    // Setting the mask to average is due to the reasoning that if the mask isn't changed, it is best that what is generated is not just a copy of one adult (hence, not setting it to Parent 1 or 2)
    for (int i = 0; i < OPTIM_VARS; i++)
    {
        // TODO: is this the best initial setting? Does this even matter?
        mask.push_back(AVG);
    }

    // Vector that determines which parents will be paired up to generate children
    // it will contain indexes of adults that are deemed to be parents
    // it will be shuffled to generate random parent pairings
    // NOTE: this will only be effective after oldAdults has been sorted (with the best Adults at the front)
    std::vector<int> parentPool;

    // Fill the parentPool index with the number of indexes desired for survivors
    for (int i = 0; i < parents.size(); i++)
    {
        parentPool.push_back(i);
    }

    // Shuffle the parentPool mapping array
    // This will generate a random list of indicies, which can then be used to map to random Adults within oldAdult's survivor range to pair parents
    // This will effectively generate random parent pairs
    std::shuffle(parentPool.begin(), parentPool.end(), rng);

    // Tracks the number of new individuals that have been generated
    // Used primarily to make sure that no children are overwritten in newChildren; indexing in generatechildrenPair, should equal childrenToGenerate when finished
    // Starts as 0 beacuse no new children have been generated
    int numNewChildren = 0;

    // Tracks the the iteration through the parentPool vector; determines which parents are passed in to the generate children function
    // Will also help track if parentPool needs to be reshuffled to generate new parent pairs
    int parentIndex = 0;

    // Toggles new pairings of parents without needing to reshuffle
    //       (parent offeset = 0: 0,1 2,3 3,4...) (parent offset = 1: 1,2 3,4, 5,6...)
    // Help keep track of how many times we have looped
    int parentOffset = 0;

    // Generate new children until we have generated enough
    // NumNewChildren will update every time generateChildrenPair is called
    // childrenToGenerate is the number of children that need to be generated
    while (numNewChildren < childrenToGenerate)
    {
        // At this point, parentIndex/+1 could possibly have crossed the size of parentPool

        // Check to see if we can calculate new pairs (else will be modifying parentIndex)
        if (parentIndex + 1 < parentPool.size())
        {
            // Generate a pair of children based on the random cross over mask method from a pair of parents
            // This will generate a set of parameters with variables randomly selected from each parent

            // Generate the base mask, with each variable being assigned to a random parent
            crossOver_wholeRandom(mask, rng);

            // Generate a pair of children based on the mask
            generateChildrenPair(parents[parentPool[parentIndex]], parents[parentPool[parentIndex + 1]], newChildren, childrenToGenerate, mask, currentAnneal, rng, numNewChildren, generation, cConstants);

            // Generate a pair of children from a pair of parents based on the average mask method
            // This mask method will generate a set of parameters based on the average of parameters from each parent

            // Set the mask to be be average
            crossOver_average(mask);
            // Generate a pair of children based on the mask
            generateChildrenPair(parents[parentPool[parentIndex]], parents[parentPool[parentIndex + 1]], newChildren, childrenToGenerate, mask, currentAnneal, rng, numNewChildren, generation, cConstants);

            // Generate a pair of children from a pair of parents based on the bundled random mask method
            // This mask method will generate parameters with variables pulled from random parents
            // This is different from the wholeRandom method in that the thrust variables (gamma, tau, and coast) is taken from the same parent

            // Generate the base mask, each variable (or set of variables) is set to a random parent
            crossOver_bundleVars(mask, rng);
            // Generate a pair of children based on the mask
            generateChildrenPair(parents[parentPool[parentIndex]], parents[parentPool[parentIndex + 1]], newChildren, childrenToGenerate, mask, currentAnneal, rng, numNewChildren, generation, cConstants);

            // Add two to parentIndex to account for the variable tracking pairs of parents, not just one parent's index
            parentIndex += 2;
        }
        else
        { // ParentIndex has crossed the boundary of parentPool size
            if (parentOffset == 1)
            {
                std::cout << "\nIn generateChildrenFromCrossover() we've cycled through the parentPool twice\n";

                // Force new pairs
                std::shuffle(parentPool.begin(), parentPool.end(), rng);

                parentOffset = 0;
                parentIndex = 0;
            }
            else
            {
                parentIndex = 1;
                // Forces pairs to be off by 1
                parentOffset = 1;
            }
        }
        //A new pair was computed or we reset parentOffset
    }   //Repeat loop
    //std::cout << "\n_-_-_-_-_-_-_-_-_-Number of children generatated: " << numNewChildren << "_-_-_-_-_-_-_-_-_-\n";
}

//Generate mutated children from duplicate adults
void generateChildrenFromMutation(std::vector<Adult> & duplicates, Child* newChildren, const int & startingIndex, std::mt19937_64 & rng, const double& currentAnneal, const int & generation, const cudaConstants* cConstants) {

    if (cConstants -> num_individuals - startingIndex > duplicates.size()) 
    {
        std::cout << "\nERROR: Scott says this is an error (beginning of mutate children),\n\tDuplicates is too small!\n";
    }

    //Go through the duplicates vector and use the adults' info to generate children into newChildren
    //NOTE: The index that is used to generate the child within newChildren depends on the size of the duplicates vector
    //          The total number of new children that need to be created between both functions is num_individuals
    //          generateChildrenFromCrossover has already filled newChildren with startingIndex number of children
    //          Thus, the starting index until num_individuals is safe to be modified
    //bool dupe = true;
    for (int i = 0; i < cConstants->num_individuals - startingIndex; i++) {
        //Assign the starting parameters of a corresponding duplicate to a child
        newChildren[startingIndex + i] = Child(duplicates[i].startParams, cConstants, generation);

        //Now that the necessary child has been assigned starting parameters, go though and heavily mutate their parameters
        //      Note: the mutation factor is set by duplicate_mutation_factor
        newChildren[startingIndex + i].startParams = mutate(newChildren[startingIndex + i].startParams, rng, currentAnneal, cConstants, generation, cConstants->mutation_amplitude, cConstants->duplicate_mutation_chance);
    }
}
