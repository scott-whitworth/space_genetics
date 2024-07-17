// Creates the next pool to be used in the optimize function in opimization.cu
void newGeneration (std::vector<Adult> & oldAdults, std::vector<Adult> & newAdults, const int & generation, std::mt19937_64 & rng, const cudaConstants* cConstants, GPUMem & gpuValues) {

    //Vector that will hold the adults who are potential parents
    //The criteria for being a parent is being in the top survivor_count number of adults in the oldAdult pool
    std::vector<Adult> parents; 

    parents.clear();

    //separates oldAdults into parents (full of best N/4 individuals) 
    fillParents(oldAdults, parents, generation, cConstants);

    //Import number of new individuals the GPU needs to fill its threads
    //Create a newChildren function to fill in with children generated from oldAdults
    Child* newChildren = new Child[cConstants->num_individuals]; 

    //uses parents to make children for the next generation
    makeChildren(parents, newChildren, generation, rng, cConstants);
    
    //Initialize variables needed for callRK
    double timeInitial = 0;
    double calcPerS = 0;

    // Each child represents an individual set of starting parameters 
    // GPU based runge kutta process determines final position and velocity based on parameters
    // Will fill in the final variables (final position & speed, posDiff, speedDiff) for each child
    //    stepSize -> (cConstants->orbitalPeriod / cConstants->max_numsteps): 
    //                good first guess as to step size, if too small/large RK will always scale to error tolerance
    callRK(calcPerS, newChildren, cConstants, gpuValues, timeInitial);

    //Now that the children have been simulated, convert the children into adults
    //First, it will calculate the right pos and speed diffs for the children
    //This will also put the converted children into the newAdults vector
    convertToAdults(newAdults, newChildren, cConstants); 

    //Free the pointer's memory
    delete[] newChildren; 
}

void fillParents(std::vector<Adult> & oldAdults, std::vector<Adult> & parents, const int & generation, const cudaConstants* cConstants){
    //Clearing parents before separating the parents from a vector ensures parents
    parents.clear();

    //If statement will check if the generation is 0 
    //      If so, it will assign all the surviving adults to parents
    //      If not, it will decide which vector based on progress
    //      This is done because, only in generation 0, does the adults in oldAdults not have a calculated progress
    if (generation == 0) {
        //Assign all surviving adults to parents
        for (int i = 0; i < cConstants->survivor_count; i++) {
            parents.push_back(oldAdults[i]);
        }   
    }else{

        //sort all the oldAdults from best to worst, so the best will be chosen as survivors/parents
        // std::sort(oldAdults.begin(), oldAdults.end(), rankRaritySort);
        mainSort(oldAdults, cConstants, oldAdults.size());
        
        //Iterate through the best of oldAdults and sort the individuals into parents
        for (int i = 0; i < cConstants->survivor_count; i++)
        {
            if (oldAdults[i].errorStatus == VALID || oldAdults[i].errorStatus == DUPLICATE)
            {
                parents.push_back(oldAdults[i]);
            }
        }

        //Check to see if the size of parents is less then survivor_count
        //This means there were adults with additional errors that made it into the survivor range
        if (parents.size() < cConstants->survivor_count)
        {
            std::cout << "\nThere are less adults in parents than desired!\n"; 
        }
    }

}

//
void makeChildren(std::vector<Adult> & parents, Child * newChildren, const int & generation, std::mt19937_64 & rng, const cudaConstants* cConstants){

    //Variable that tracks the number of children that needs to be generated via the crossover method
    //Set to the number of children that needs to be generated (num_individuals)
    int childrenFromCrossover = cConstants->num_individuals;

    //Generate children with crossovers using parents
    generateChildrenFromCrossover(parents, newChildren, childrenFromCrossover, rng, generation, cConstants);
    

}

//Function that will convert the generated children into adults
//Will also transfer the created adults into the newAdult vector
void convertToAdults(std::vector<Adult> & newAdults, Child* newChildren, const cudaConstants* cConstants) {
    //Clear newAdults as a percaution against newAdults becoming too big
    newAdults.clear();

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
        }//if it is not a nan, the status has already been made valid or sun_error in rk4CUDASim

        //Now that the status has been determined, there is enough information to set pos and speed diffs
        //The two functions will look at the child's errorStatus and set the diffs based on that
        newChildren[i].getPosDiff(cConstants);
        newChildren[i].getSpeedDiff(cConstants);
        newChildren[i].getHorzVelDiff(cConstants);
        newChildren[i].getVertVelDiff(cConstants);
        newChildren[i].getOrbitPosDiff(cConstants);
        newChildren[i].getOrbitSpeedDiff(cConstants);

        for (int j = 0; j < cConstants->missionObjectives.size(); j++) {
            newChildren[i].objTargetDiffs[j] = abs(newChildren[i].getParameters(cConstants->missionObjectives[j]) - cConstants->missionObjectives[j].target);
        }

        //Get progress now that all of the outputs are ready
        newChildren[i].getProgress(cConstants);
        
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
void createFirstGeneration(std::vector<Adult>& oldAdults, const cudaConstants* cConstants, std::mt19937_64 rng, const int & generation, GPUMem & gpuValues, ReferencePoints & refPoints){
    //int generation = 0; //these are loaded from a file or randomly, so set gen to 0
    Child* initialChildren = new Child[cConstants->num_individuals]; 

    //Try to open the final allAdults excel file from the last run
    //  Note: the time seeds are subtrracted by 100 because the main difference between runs is the time seed is incremented by 100, so this line will get the previous run's file
    std::string prevIndPath = "../Output_Files/" + std::to_string(static_cast<int>(cConstants->time_seed-100)) + "_1/" + std::to_string(static_cast<int>(cConstants->time_seed-100)) + "-AllAdults-gen#" + std::to_string(cConstants->max_generations) + ".csv";
    std::fstream prevIndividuals(prevIndPath);

    std::cout << "\nChecked Path: " << prevIndPath << "\n";
    //If this does not open, there are two likely scenerios
    //  1) This is the first run, so the folder for the previous run doesn't exist
    //  2) The previous run converged, so the last allAdults excel file isn't the maxGeneration one
    //For both of these scenerios, we need to randomly generate initial parameters

    // See if there should be a random start
    //  If carryover_individuls is 0, it means that every run should have a random start
    //  If the fream isn't good, it means the file didn't open properly (see above for normal reasons why) so the initial values should be randomly generated
    if (cConstants->carryover_individuals == 0 || !prevIndividuals.good()) {
        std::cout << "\nStarting with random parameters!\n";
        // individuals set to randomly generated, but reasonable, parameters
        for (int i = 0; i < cConstants->num_individuals; i++) { 
            initialChildren[i] = Child(randomParameters(rng, cConstants), cConstants, generation, 0);
        }
    }
    // Get the initial conditions from the previous run
    else {   
        std::cout << "\nStarting with previous run's parameters!\n";
        // sort the data into 2 dimensions
        // each row is an individual
        // each column is one of the starting parameters
        std::string startParamVal;
        // arrayCPU needs to be updated to handle the fact that OPTIM_VARS may be flexible
        // double arrayCPU[OPTIM_VARS][cConstants->carryover_individuals];
        std::vector<std::vector<double>> baseParams (OPTIM_VARS, std::vector<double>(cConstants->carryover_individuals, 0));

        //String that will store the string of the line
        std::string adultInfo;
        
        //Create pivots to parse through line info from allAdults
        int startPivot, endPivot;

        //Initial getline to clear the header row
        std::getline(prevIndividuals, adultInfo);

        //Tracker will count how many columns need to be skipped before encounering the input variables
        int parseCount = 0;

        //Set the initial start pivot to the initial character
        startPivot = 0;
        //Set the intial end pivot to be the first comma after the start pivot
        endPivot = adultInfo.find(",", startPivot);

        //To move through the columns to find where the input variables start
        //  First check to see if the substring between the pivots is 'gamma0' (the first input variable printed)
        while ((adultInfo.substr(startPivot, endPivot-startPivot)) != "gamma0") {
            //the column is not the start of the input vars, select the next pivot points
            //Set the start pivot as the next character after the end pivot
            startPivot = endPivot + 1;
            //Set the next end pivot as the next comma
            endPivot = adultInfo.find(",", startPivot);

            //Add one to parse count to signify a column being skipped
            parseCount++;
        }

        //subtract one from parse count to account for 0 indexing
        parseCount--;

        //for loop will continue until the number of parameter sets wanted have been pulled from prevIndividuals
        for (int i = 0; i < cConstants->carryover_individuals; i++) {
            //Get the new line and check to see if the end of the csv has been reached
            if ( std::getline(prevIndividuals, adultInfo) ) {
                //std::cout << "\nPulled line:\n\t" << adultInfo << "\n";

                //A new adult has been pulled, get its starting params
                //Reset pivots
                startPivot = endPivot = 0;

                //Set the initial start pivot to the character after the first comma
                startPivot = adultInfo.find(",") + 1;

                //For loop will get all of the paramters for the line
                //Only need to pase for the number of skips plus the number of input parameters that need to be pulled
                for (int j = 0; j < (OPTIM_VARS + parseCount); j++) {
                    
                    //Get the next end pivot
                    //it is the first comma after the start pivot
                    endPivot = adultInfo.find(",", startPivot);

                    //Check to see if this is a start parameter
                    //  Checked by seeing if the current column is past the ones that need to be skipped
                    //  The check uses j+1 to account for 0 index
                    if (j >= parseCount) {
                        //The input param is the substring between the two pivots 
                        startParamVal = adultInfo.substr(startPivot, endPivot-startPivot);

                        //Push the value to the arrayCPU array
                        //The order of the params in the csv file is the same order as the OPTIM_VARS array, so no need to parse
                        baseParams[j-parseCount][i] = std::stod(startParamVal);
                    }

                    //Resetting the start pivot so the next parameter can be pulled
                    startPivot = endPivot + 1;
                }

            }
            else {
                //The csv ended before the requested number of starting params were pulled
                //This is an error, since carryover_individuals should never be larger than num_individuals
                //  num_individuals itself would normally only half as large as the number of adults in allAdults
                //  so something went very wrong (likely many errors at the end of the previous run) if the code is here
                std::cout << "\nERROR: too few rows in previous run's final allAdults file!\n";
                break;
            }
        }

        //Counter to cycle through the saved parameters
        int paramNum = 0;
        
        //While loop that will create the children from the saved parameters
        for (int i = 0; i < cConstants->num_individuals; i++) {
            //Get the stored parameters for this individual
            double tripTime = baseParams[TRIPTIME_OFFSET][paramNum];
            double alpha = baseParams[ALPHA_OFFSET][paramNum];
            double beta = baseParams[BETA_OFFSET][paramNum];
            double zeta = baseParams[ZETA_OFFSET][paramNum];

            coefficients<double> testcoeff;
            for (int j = 0; j < testcoeff.gammaSize; j++) {
                testcoeff.gamma[j] = baseParams[j + GAMMA_OFFSET][paramNum];
            }

            for (int j = 0; j < testcoeff.tauSize; j++) {
                testcoeff.tau[j] =  baseParams[j + TAU_OFFSET][paramNum];
            }

            for (int j = 0; j < testcoeff.coastSize; j++) {
                testcoeff.coast[j] = baseParams[j + COAST_OFFSET][paramNum];
            }

            //Create the new parameters
            rkParameters<double> example(tripTime, alpha, beta, zeta, testcoeff);

            //Mutate the parameters
            example = mutate(example, rng, cConstants->anneal_initial, cConstants, cConstants->mutation_amplitude, cConstants->default_mutation_chance); 

            //Create a new child with the mutated parameters
            initialChildren[i] = Child(example, cConstants, generation, 0); 

            //Add one to paramNum to get the next set of parameters
            paramNum++;

            //Check to make sure paramNum doesn't go above the size of the array (checking the size of the gamma_offset column is arbitrary)
            if (paramNum >= baseParams[GAMMA_OFFSET].size()) {
                //Reset paramNum
                paramNum = 0;
            }    
        }
    }
    firstGeneration(initialChildren, oldAdults, cConstants, rng, gpuValues, refPoints);
    delete[] initialChildren;
}

//Will create the first generation of adults from random parameters so that future generations have a pre-existing base of adults
void firstGeneration(Child* initialChildren, std::vector<Adult>& oldAdults, const cudaConstants* cConstants, std::mt19937_64 rng, GPUMem & gpuValues, ReferencePoints & refPoints){
    double timeIntial = 0;
    double calcPerS  = 0;

    // Each child represents an individual set of starting parameters
    // GPU based runge kutta process determines final position and velocity based on parameters
    //Will fill in the final variables (final position & speed, posDiff, speedDiff) for each child
    callRK(calcPerS, initialChildren, cConstants, gpuValues, timeIntial);
    
    //Now that the children have been simulated, convert the children into adults
    //This will calc the objective outputs for each child
    //It will also put the converted children into the newAdults vector
    convertToAdults(oldAdults, initialChildren, cConstants);
}

//fills oldAdults with the best adults from this generation and the previous generation so that the best parents can be selected
void preparePotentialParents(std::vector<Adult>& allAdults, std::vector<Adult>& newAdults, std::vector<Adult>& oldAdults, int& numErrors, int& duplicateNum, const cudaConstants* cConstants, std::mt19937_64 rng, ReferencePoints & refPoints, const int & generation, int & marsErrors){
    //ensures this vector is empty and ready for new inputs
    std::vector<Adult>().swap(allAdults); 

    //Reset duplicate and error count variables
    numErrors = 0;
    marsErrors = 0;
    duplicateNum = 0;

    //Iterate through allAdults and find any duplicate adults
    findDuplicates(newAdults, oldAdults, cConstants);

    //get rid of any invalid or old adults, as well as adults with bad posDiffs or speedDiffs that are too small
    //check the errorStatus of all the newAdults and add them to allAdults
    for (int i = 0; i < newAdults.size(); i++){ 
        //copies all the elements of newAdults into allAdults
        addToAllAdults(newAdults, allAdults, i, numErrors, duplicateNum, marsErrors);
    }
    //check the errorStatus of all the oldAdults and add them to allAdults
    for (int i = 0; i < oldAdults.size(); i++){ //copies over all the elements of oldAdults into allAdults
        addToAllAdults(oldAdults, allAdults, i, numErrors, duplicateNum, marsErrors);
    }

    if(allAdults.size() < cConstants->num_individuals){
        std::cout << "The size of allAdults is smaller than N (num_individuals)" << std::endl;
    }else if(allAdults.size() < cConstants->survivor_count){
        std::cout << "ERROR: The size of allAdults is smaller than survivor count!" << std::endl;
    }

    //Calculate rank of each adult based on domination sort
    //* Ignore any nans at the end of allAdults
    //must be called after checking for nans and before giveDistance
    giveRank(allAdults, cConstants); //gives a rank to each adult

    //Check to see which algorithm to use
    if (cConstants->algorithm == RANK_RARITY) {
        //The user wants to use the rank-rarity method, call the necessary functions for it

        //Calculate the rarity for all Adults
        giveRarity(cConstants, rng, refPoints, allAdults);

        //Adjust the rarities to assigned the reserved scores to the best individuals
        assignReservedRarity(cConstants, allAdults);

        /* The code below is used for calculating the rarity within each rank
        //Sort the adults by rank
        std::sort(allAdults.begin(), allAdults.end(), rankSort);

        //Assign rarity for each rank
        //Create a vector which stores the adults in a specific rank
        std::vector<Adult> rankAdults(0);

        //ints to track the indexes of the first and last adults in a rank
        int start = 0, end = 0;

        //Calculate rarity for each rank
        //  Calculate rarity for the total number of ranks, which is easily found by getting the rank of the last inidividual
        while (end < allAdults.size()) {
            //Interate the end index tracker until the next rank is reached or it reaches the end
            while (allAdults[end].rank == allAdults[start].rank && end < allAdults.size()) {
                //Add the rank adult to the rankAdults vector
                rankAdults.push_back(allAdults[end]);

                //Move to the next adult
                end++;
            }
            
            //Calculate the rarity for the rank
            giveRarity(cConstants, rng, refPoints, rankAdults);

            //Put the ranked adults back into allAdults
            for (int k = 0; k < rankAdults.size(); k++) {
                allAdults[start + k] = rankAdults[k];
            }

            //Set the start to be the new end
            start = end;

            //clear the rankAdults vector for the next rank
            rankAdults.clear();
        }
        */
    }
    else {
        //The user wants the rank-distance method or the algorithm is unspecified (use rank-distance by default)
        giveDistance(allAdults, cConstants); //gives a distance to each adult
    }

    //sorts allAdults using preferred sort
    mainSort(allAdults, cConstants, allAdults.size());

    oldAdults.clear(); //empties oldAdults so new values can be put in it
    int counter = 0; //a variable that ensures we put the correct number of adults in oldAdults

    //copies the best adults from allAdults into oldAdults (should be half of allAdults that are copied over)
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
}

void addToAllAdults(std::vector<Adult> & adultPool, std::vector<Adult> & allAdults, const int & index, int& numErrors, int& duplicateNum, int& marsErrors){
    //if we ever want to eliminate duplicates again, this is the place to do it
    
    if(adultPool[index].errorStatus != VALID && adultPool[index].errorStatus != DUPLICATE) { 
        //tallies the number of nans in allAdults by checking if the adult being passed into newAdult is a Nan or not
        numErrors++;
        
        if (adultPool[index].errorStatus == MARS_ERROR){
            marsErrors++;
        }

    }else if(adultPool[index].errorStatus == DUPLICATE){
        duplicateNum++;
        // allAdults.push_back(adultPool[index]);//remove this if we want to remove duplicates
    }
    else{
        allAdults.push_back(adultPool[index]);
    }
}