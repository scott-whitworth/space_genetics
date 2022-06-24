// Creates the next pool to be used in the optimize function in opimization.cu
void newGeneration (std::vector<Adult> & oldAdults, std::vector<Adult> & newAdults, const double & annealing, const int & generation, std::mt19937_64 & rng, const cudaConstants* cConstants) {

    //Vector that will hold the adults who are potential parents
    //The criteria for being a parent is being in the top survivor_count number of adults in the oldAdult pool and not being a duplicate (distance > 0)
    //      Duplicate adults will generate children in a separate fashion
    std::vector<Adult> parents; 

    //Vector for duplicates, based on the same criteria as above
    std::vector<Adult> duplicates;

    parents.clear();
    duplicates.clear();

    //separates oldAdults into parents (full of unique individuals) and duplicates (anything whose posDiff and speedDiff are about the same as an Adult in parents)
    separateDuplicates(oldAdults, parents, duplicates, generation, cConstants);

    //Import number of new individuals the GPU needs to fill its threads
    //Create a newChildren function to fill in with children generated from oldAdults
    Child* newChildren = new Child[cConstants->num_individuals]; 

    //uses parents and duplicates to make children for the next generation
    makeChildren(parents, duplicates, newChildren, annealing, generation, rng, cConstants);
    
    //Initialize variables needed for callRK
    //TODO: there is likely a better solution for this
    double timeInitial = 0;
    double calcPerS = 0;

    // Each child represents an individual set of starting parameters 
    // GPU based runge kutta process determines final position and velocity based on parameters
    //Will fill in the final variables (final position & speed, posDiff, speedDiff) for each child
    //    stepSize -> (cConstants->orbitalPeriod / cConstants->max_numsteps): 
    //                good first guess as to step size, if too small/large RK will always scale to error tolerance
    //TODO: Perhaps we should just import cCOnstants and newChildren into callRk, since most of the arguments come from cConstants regardless
    callRK(cConstants->num_individuals, cConstants->thread_block_size, newChildren, timeInitial, (cConstants->orbitalPeriod / cConstants->max_numsteps), cConstants->rk_tol, calcPerS, cConstants); 

    //Now that the children have been simulated, convert the children into adults
    //First, it will calculate the right pos and speed diffs for the children
    //This will also put the converted children into the newAdults vector
    convertToAdults(newAdults, newChildren, cConstants); 

    //Free the pointer's memory
    delete[] newChildren; 
}

//separates the duplicates from the unique individuals (parents)
void separateDuplicates(std::vector<Adult> & oldAdults, std::vector<Adult> & parents, std::vector<Adult> & duplicates, const int & generation, const cudaConstants* cConstants){
    //If statement will check if the generation is 0 
    //      If so, it will assign all the surviving adults to parents
    //      If not, it will decide which vector based on cost
    //      This is done because, only in generation 0, does the adults in oldAdults not have a calculated cost
    if (generation == 0) {
        //Assign all surviving adults to parents
        for (int i = 0; i < cConstants->survivor_count; i++) {
            parents.push_back(oldAdults[i]);
        }   
    }else{
        //TODO: This should be its own function (which should make unit testing easier) - Done?
        //sort all the oldAdults by rankDistance, so the best will be chosen as survivors/parents
        std::sort(oldAdults.begin(), oldAdults.end(), rankDistanceSort);
        
        //Iterate through the best of oldAdults and sort the individuals into either the parent or the duplicate vectors
        for (int i = 0; i < cConstants->survivor_count; i++)
        {
            if (oldAdults[i].errorStatus == VALID)
            {
                parents.push_back(oldAdults[i]);
            }
            else if (oldAdults[i].errorStatus == DUPLICATE)
            {
                duplicates.push_back(oldAdults[i]);
            }
        }

        //Check to see if the size of parents + the size of duplicates is less then survivor_count
        //This means there were adults with additional errors that made it into the survivor range
        if (parents.size() + duplicates.size() < cConstants->survivor_count)
        {
            std::cout << "\nThere are less adults in parents & duplicates than desired!\n"; 
        }
    }

}

//
void makeChildren(std::vector<Adult> & parents, std::vector<Adult> & duplicates, Child * newChildren, const double & annealing, const int & generation, std::mt19937_64 & rng, const cudaConstants* cConstants){
    //TODO: start of next function - DONE
    //Variable that tracks the number of children that needs to be generated via the crossover method
    //Set to the number of children that needs to be generated (num_individuals) minus the number of children that will be generated from duplicates via heavy mutation (size of duplicates)
    int childrenFromCrossover = cConstants->num_individuals - duplicates.size();    

    //Generate children with crossovers using non-duplicate parents
    generateChildrenFromCrossover(parents, newChildren, childrenFromCrossover, rng, annealing, generation, cConstants);

    //See if there are any duplicate adults before generating children from mutations
    if(duplicates.size() > 0){
        //Generate the rest of the children using heavy mutation of duplicate adults
        generateChildrenFromMutation(duplicates, newChildren, childrenFromCrossover, rng, annealing, generation, cConstants); 
    }
}

//Function that will convert the generated children into adults
//Will also transfer the created adults into the newAdult vector
void convertToAdults(std::vector<Adult> & newAdults, Child* newChildren, const cudaConstants* cConstants) {
    //TODO: Add a clear to newAdults, just in case
    //      This should all be documented clearly in the header

    //Iterate through the simulated children to determine their error status (if it hasn't been set by callRK already) and calculate their pos and speed differences
    for (int i = 0; i < cConstants->num_individuals; i++)
    {
        //Check to see if nans are generated in the finalPos elements
        if (  isnan(newChildren[i].finalPos.r) ||
              isnan(newChildren[i].finalPos.theta) ||
              isnan(newChildren[i].finalPos.z) ||
              isnan(newChildren[i].finalPos.vr) ||
              isnan(newChildren[i].finalPos.vtheta) ||
              isnan(newChildren[i].finalPos.vz)  ) {
            //Mark the child with the nan variables with the nan_error flag
            newChildren[i].errorStatus = NAN_ERROR;
        }//if it is not a nan, the status has already been made valid or sun_error in rk4SimpleCuda

        //Now that the status has been determined, there is enough information to set pos and speed diffs
        //The two functions will look at the child's errorStatus and set the diffs based on that
        newChildren[i].getPosDiff(cConstants);
        newChildren[i].getSpeedDiff(cConstants); 
    }
    
    //Iterate through the newChildren vector to add them to the newAdult vector
    //iterates until num_individuals as that is the number of new children needed to generate
    //      This accounts for the potential that the number of children in newChildren is slightly larger than num_individuals
    for (int i = 0; i < cConstants->num_individuals; i++)
    {
        //Fill in the newAdults vector with a new adult generated with a completed child
        newAdults.push_back(Adult(newChildren[i]));
    }
}


//Creating a generation either randomly or based on values from a file
void createFirstGeneration(std::vector<Adult>& oldAdults, const cudaConstants* cConstants, std::mt19937_64 rng, const int & generation){
    //int generation = 0; //these are loaded from a file or randomly, so set gen to 0
    Child* initialChildren = new Child[cConstants->num_individuals]; 
    // Initialize individuals randomly or from a file
    if (cConstants->random_start) {
        // individuals set to randomly generated, but reasonable, parameters
        for (int i = 0; i < cConstants->num_individuals; i++) { 
            initialChildren[i] = Child(randomParameters(rng, cConstants), cConstants, generation);
        }
    }
    // Read from file using cConstants initial_start_file_address to get path
    else {
        // **Might be depreciated, not tested summer 2020**
        // Sets inputParameters to hold initial individuals based from file optimizedVector.bin
        const int numStarts = 14; // the number of different sets of starting parameters in the input file
        std::ifstream starts;
        starts.open(cConstants->initial_start_file_address, std::ifstream::in|std::ios::binary); // a file containing the final parameters of converged results from CPU calculations        

        // sort the data into 2 dimensions
        // one row is one set of starting parameters
        // each column is a specific variable:
        double startDoubles;
        // arrayCPU needs to be updated to handle the fact that OPTIM_VARS may be flexible
        double arrayCPU[numStarts][OPTIM_VARS];
        
        for (int i = 0; i < OPTIM_VARS; i++) { // rows
            for (int j = 0; j < numStarts; j++) { // columns
                starts.read( reinterpret_cast<char*>( &startDoubles ), sizeof startDoubles );
                arrayCPU[j][i] = startDoubles;
            }
        }
        starts.close();

         // set every thread's input parameters to a set of final values from CPU calculations for use as a good starting point
        for (int i = 0; i < cConstants->num_individuals; i++) {
            int row = rng() % numStarts; // Choose a random row to get the parameters from

            double tripTime = arrayCPU[row][TRIPTIME_OFFSET];
            double alpha = arrayCPU[row][ALPHA_OFFSET];
            double beta = arrayCPU[row][BETA_OFFSET];
            double zeta = arrayCPU[row][ZETA_OFFSET];

            coefficients<double> testcoeff;
            for (int j = 0; j < testcoeff.gammaSize; j++) {
                testcoeff.gamma[j] = arrayCPU[row][j + GAMMA_OFFSET];
            }

            for (int j = 0; j < testcoeff.tauSize; j++) {
                testcoeff.tau[j] =  arrayCPU[row][j + TAU_OFFSET];
            }

            for (int j = 0; j < testcoeff.coastSize; j++) {
                testcoeff.coast[j] = arrayCPU[row][j + COAST_OFFSET];
            }

            rkParameters<double> example(tripTime, alpha, beta, zeta, testcoeff); 

            initialChildren[i] = Child(example, cConstants, generation);
        }
    }
    firstGeneration(initialChildren, oldAdults, cConstants);
    delete[] initialChildren;
}

//Will create the first generation of adults from random parameters so that future generations have a pre-existing base of adults
void firstGeneration(Child* initialChildren, std::vector<Adult>& oldAdults, const cudaConstants* cConstants){
    double timeIntial = 0;
    double calcPerS  = 0;

     // Each child represents an individual set of starting parameters 
        // GPU based runge kutta process determines final position and velocity based on parameters
    //Will fill in the final variables (final position & speed, posDiff, speedDiff) for each child
    //TODO: Perhaps we should just import cCOnstants and newChildren into callRk, since most of the arguments come from cConstants regardless
    callRK(cConstants->num_individuals, cConstants->thread_block_size, initialChildren, timeIntial, (cConstants->orbitalPeriod / cConstants->max_numsteps), cConstants->rk_tol, calcPerS, cConstants); 

    //Now that the children have been simulated, convert the children into adults
    //This will calc the posDiff and speedDiff for each child
    //It will also put the converted children into the newAdults vector
    convertToAdults(oldAdults, initialChildren, cConstants); 
}

//fills oldAdults with the best adults from this generation and the previous generation so that the best parents can be selected
void preparePotentialParents(std::vector<Adult>& allAdults, std::vector<Adult>& newAdults, std::vector<Adult>& oldAdults, int& numNans, int& duplicateNum, const cudaConstants* cConstants, const int & generation, const double& currentAnneal){
    std::vector<Adult>().swap(allAdults); //ensures this vector is empty and ready for new inputs
    numNans = 0;
    duplicateNum = 0;

    //TODO: Making sure there are no errors within allAdults is temporary,
    //       eventually we need to make sure they are at the end of the rankDistanceSort but still include them within allAdults

    //TODO: We should figure out if we want the size of allAdults (or the other vectors) to not be their expected size
    //      The following code may result in allAdults being less than 2*n (or 2*num_individuals)
    //      There are good reasons either way for this choice
    // 
    //      We should for sure be checking to make sure allAdults never gets below survivor size / num_individuals
    //      

    //countStatuses(newAdults, generation);
    //countStatuses(oldAdults, generation);
    //Iterate through allAdults and find any duplicate adults
    findDuplicates(newAdults, oldAdults, cConstants, currentAnneal);

    //get rid of any invalid or old adults, as well as adults with bad posDiffs or speedDiffs that are too small
    eliminateBadAdults(allAdults, newAdults, oldAdults, numNans, duplicateNum, cConstants, generation);

    

    //countStatuses(allAdults, generation);

    //Calculate rank of each adult based on domination sort
    //* Ignore any nans at the end of allAdults
    //must be called after checking for nans and before giveDistance
    //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE GIVE RANK-_-_-_-_-_-_-_-_-_\n\n";
    giveRank(allAdults, cConstants); //gives a rank to each adult
    //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE GIVE DISTANCE-_-_-_-_-_-_-_-_-_\n\n";
    giveDistance(allAdults, cConstants); //gives a distance to each adult
    //std::cout << "\n\n_-_-_-_-_-_-_-_-_-TEST: PRE GIVE SORT-_-_-_-_-_-_-_-_-_\n\n";
    std::sort(allAdults.begin(), allAdults.end(), rankDistanceSort); //sorts allAdults using rankDistance sort

    oldAdults.clear(); //empties oldAdults so new values can be put in it
    int counter = 0; //a variable that ensures we put the correct number of adults in oldAdults
    //copies the best adults from allAdults into oldAdults (should be half of allAdults that are copied over)
    //TODO:: Make sure dividing num_individuals / 2 is correct or if we should divide by something else

    while (counter < cConstants->num_individuals && counter < allAdults.size()){
        oldAdults.push_back(allAdults[counter]);
        counter++;
    }

    //error message prints if somehow there are not enough adults to fill oldIndividuals to the size it should be  
    if(counter == allAdults.size() && counter < cConstants->num_individuals){ //TODO: May or may not want this. If used before oldAdults and newAdults were both filled once, delete error message
        std::cout << "There are not enough adults to fill oldIndividuals properly" << std::endl;
    }

    //Clear the newAdult vector so that it is ready to receive new children
    newAdults.clear(); 

    //countStatuses(oldAdults, generation);
}

//eliminates unwanted adults from allAdults
void eliminateBadAdults(std::vector<Adult>& allAdults, std::vector<Adult>& newAdults, std::vector<Adult>& oldAdults, int& numNans, int& duplicateNum, const cudaConstants* cConstants, const int & generation){
    double posTolerance = cConstants->pos_threshold/100;
    double speedTolerance = cConstants->speed_threshold/100;
    int ageTolerance = 500;

    for (int i = 0; i < newAdults.size(); i++){ //copies all the elements of newAdults into allAdults
        if(newAdults[i].errorStatus != VALID && newAdults[i].errorStatus != DUPLICATE) { //tallies the number of nans in allAdults by checking if the adult being passed into newAdult is a Nan or not
            numNans++;
        }else if((newAdults[i].posDiff < posTolerance && newAdults[i].speedDiff > speedTolerance) || (newAdults[i].posDiff > posTolerance && newAdults[i].speedDiff < speedTolerance)){//check if the adult is too good in one aspect
            //do nothing and don't add them to all adults
        }else if(newAdults[i].errorStatus == DUPLICATE){
            duplicateNum++;
        }
        else {
            allAdults.push_back(newAdults[i]);
        }
    }
    for (int i = 0; i < oldAdults.size(); i++){ //copies over all the elements of oldAdults into allAdults
        if(oldAdults[i].errorStatus != VALID && oldAdults[i].errorStatus != DUPLICATE){//tallies the number of nans in allAdults by checking if the adult being passed into newAdult is a Nan or not
            numNans++;
            //TODO: We should cout an error here, there should not ever be oldAdults with error status (but what about duplicates?)
        }else if(generation - oldAdults[i].birthday > ageTolerance){
            //make sure there are no adults that are too old (older than 50 generations)
            //do nothing and don't add them 
        }else if(oldAdults[i].errorStatus == DUPLICATE){
            duplicateNum++;
        }
        else {
            allAdults.push_back(oldAdults[i]);
        }
    }
    //just for debugging to make syre we are not getting rid of too many adults
    if(allAdults.size() < cConstants->num_individuals){
        std::cout << "The size of allAdults is smaller than N" << std::endl;
    }else if(allAdults.size() < cConstants->survivor_count){
        std::cout << "ERROR: The size of allAdults is smaller than N/4" << std::endl;
    }
}