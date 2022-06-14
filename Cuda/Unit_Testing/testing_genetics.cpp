bool runGeneticsUnitTests(){
// SETTING UP CUDA CONSTANTS TO BE USED BY OTHER FUNCTIONS 
    //making cudaConstants to control what is going on in genetic algorithm while still using the original functions
    cudaConstants* utcConstants = new cudaConstants(); 

    // Seed used for randomization rng things 
    // Says to set the time_seed to NONE so time(0), but that does not work so set it to 0
    utcConstants->time_seed = 0; 

    //values taken directly from the genetic config file
    utcConstants->anneal_initial = 0.10; // initial value for annealing, meant to replace the previously used calculation involving ANNEAL_MIN and ANNEAL_MAX with something more simple
    utcConstants->anneal_final = 1.0e-7;   // final value for annealing, anneal cannot be reduced beyond this point
    utcConstants->anneal_factor = 0.75;  // factor by which annealing is multiplied with when there is no change in the best individual over 100 generations

    // The percentage for probability of mutating a gene in a new individual, called iteratively to mutate more genes until the check fails
    // Starting this off at 0 because to ensure that there are no mutations initially to verify that children are generated as expected
    // Changes to 1.0 later in the code so that we can ensure mutation is working
    utcConstants->mutation_rate = 0; 
   
    // Used in mutate(), affects the scale of change for the respective parameter values, in conjunction with annealing
    // Represents max bounds of mutation, mutation will never be +/- this value
    utcConstants->gamma_mutate_scale = 3.14159; 
    utcConstants->tau_mutate_scale = 1.570795; 
    utcConstants->coast_mutate_scale = 3.14159;
    utcConstants->triptime_mutate_scale = 1.0* SECONDS_IN_YEAR;
    utcConstants->zeta_mutate_scale = 1.570795;
    utcConstants->beta_mutate_scale = 1.570795;
    utcConstants->alpha_mutate_scale = 3.14159;

    // 0 is for no thruster, 1 is for NEXT ion thruster
    // starts with no thruster to test that, then will be changed to have a thruster later on when necessary
    utcConstants->thruster_type = 0; 

    // Number of individuals in the pool -> chose to make this 10 so a generation size is managable
    // Additionally a size of 10 shows what happens if a the number of children generated is not a factor of the needed generation size
    // (E.g. 8 does not evenly divide 10 -> two pairings of parents will each produce 8 kids, which is 16 total and 16 > 10)
    utcConstants->num_individuals = 10; 

    // Number of survivors selected, every pair of survivors creates 8 new individuals 
    // Chose 3 because that is approximately a quarter of 10
    // This also allows us to see if the shuffling and new parent pairing works as expected
    // This will likely be changed in other sections of the code as well
    utcConstants->survivor_count = 3;  

    //Flag for what type of landing you want. Current modes: "soft"=1, "hard"=2
    //testing for soft because as of 2022 that is the main mission type we are exploring 
    //additionally, mission type should have little impact on the genetics algorithms
    utcConstants->missionType = 1; 

    //ALL THE STUFF SO EARTH'S INITIAL POSITION CAN BE CALCULATED
    utcConstants->triptime_max=2.0* SECONDS_IN_YEAR;
    utcConstants->triptime_min=1.0* SECONDS_IN_YEAR;
    utcConstants->timeRes=3600; // Earth Calculations Time Resolution Value
    utcConstants->v_escape = 2162.4/AU; //sqrt of DART Mission c3energy
    //Bennu config values
    utcConstants->r_fin_earth=9.857045197029908E-01;
    utcConstants->theta_fin_earth=1.242975503287042;
    utcConstants->z_fin_earth=-4.332189909686674E-05;
    utcConstants->vr_fin_earth=-1.755004992027024E-09;
    utcConstants->vtheta_fin_earth=2.019592304815492E-07;
    utcConstants->vz_fin_earth=-2.312712519131594E-12;
    // Various values that impact runge kutta
    utcConstants->rk_tol=1e-12;
    utcConstants->doublePrecThresh=1e-12;
    utcConstants->GuessMaxPossibleSteps=1000000; //TODO: Long term, I think we need to get rid of this, we use max_numsteps instead
    utcConstants->max_numsteps=2500;
    utcConstants->min_numsteps=400;

// SETTING A RANDOM NUMBER GENERATOR (rng) TO BE USED BY FUNCTIONS
    // This rng object is used for generating all random numbers in the genetic algorithm, passed in to functions that need it
    std::mt19937_64 rng(utcConstants->time_seed);

    launchCon = new EarthInfo(utcConstants); //VERY complicated part of the code with some possibility of errors -> just needed for 

// CALLING THE DIFFERENT UNIT TESTING ALGORITHMS
    bool allWorking = true;
    if (firstParentsTest(utcConstants)){
        cout << "PASSED: Children can be converted to adults that can be sorted" << endl;
    }
    else{
        cout << "FAILED: Children cannot be converted to adults that can be sorted" << endl;
        allWorking = false;
    }
    if (createMasks(rng, true)){
        cout << "PASSED: All three masks were successfully generated as expected" << endl;
    }
    else{
        cout << "FAILED: Not all the masks were successfully generated" << endl;
        allWorking = false;
    }
    if (makeChildrenWithDifferentMethods(rng, utcConstants)){
        cout << "PASSED: Successfully made children using the three different methods" << endl;
    }
    else{
        cout << "FAILED: Could not successfully make children using the three different methods" << endl;
        allWorking = false;
    }

    delete[] utcConstants;
    delete launchCon;
    return allWorking;
}

//returns true if the first generation is generated
bool firstParentsTest(const cudaConstants * utcConstants){
    std::vector<Adult> parents;
    Child* theChildren = new Child[utcConstants->num_individuals];

    //Made 10 children and set them with semi-random position difference and speed difference values that can easily be determined are either true or false
    //These numbers are semi-random multiples of 10 and were entered in no specific order at this point
    //This sort of simulates when the initial generation of oldAdults are created using random parameters 
    //(but these don't change from run to run and are easier to do compare)
    theChildren[0] = Child(150, 20);
    theChildren[1] = Child (120, 90);
    theChildren[2] = Child(180, 30);
    theChildren[3] = Child (20, 70);
    theChildren[4] = Child(110, 40);
    theChildren[5] = Child (30, 90);
    theChildren[6] = Child(50, 120);
    theChildren[7] = Child (220, 970);
    theChildren[8] = Child(20, 120);
    theChildren[9] = Child (340, 90);

    //Ranks based on my net dominations (times it dominates others - times it is dominated by others) -> dominations based on output of dominationCheckTest
    //Below is a list of the parameters of the above children and which rank they should belong to based on their position difference and speed differences
    //Rank 1: (20, 70) [NET DOMINATIONS: 6]
    //Rank 2: (150, 20); (110, 40); (30, 90) [NET DOMINATIONS: 3]
    //Rank 3: (180, 30); (20, 120) [NET DOMINATIONS: 1]
    //Rank 4: (120, 90) [NET DOMINATIONS: -1]
    //Rank 5: (50, 120) [NET DOMINATIONS: -2]
    //Rank 6: (340, 90) [NET DOMINATIONS: -6]
    //Rank 7: (220, 970) [NET DOMINATIONS-7]

    convertToAdults(parents, theChildren, utcConstants);
    wrongWayToRank(parents);
    for (int i = 0; i < parents.size(); i++){
        cout << "(" << parents[i].posDiff << "," << parents[i].speedDiff <<"): " << parents[i].unitTestingRankDistanceStatusPrint() << " / ";
    }
    cout << endl;

    sortaGiveDistance(parents);
    std::sort(parents.begin(), parents.end(), rankDistanceSort);

    //not exactly sure what the output should be here... I'm not entirely sure how giveDistance calculates precise distances to give things
    for (int i = 0; i < parents.size(); i++){
        cout << "(" << parents[i].posDiff << "," << parents[i].speedDiff <<"): " <<parents[i].unitTestingRankDistanceStatusPrint() << " / ";
    }
    cout << endl;

    delete[] theChildren;

    if (parents.size() == genSize){
        return true;
    }
    else{
        cout << "Could not generate turn children into parents to create a parent generation" << endl;
        return false;
    }

}


void wrongWayToRank(std::vector<Adult> & newAdults){
    std::vector<std::vector<int>> arrayOfDominations(genSize);  //a 2D vector that is genSize long
    for (int i = 0; i < genSize; i++){ //if there are any inputs in it so far, gives an error message
        if (arrayOfDominations[i].size() != 0){
            cout << "ERROR" << endl;
        }
    }
    int thisWasDominatedXTimes[genSize]; //an array holding how many times any singular index was dominated
    for (int h = 0; h < genSize; h++){ //sets the initial number of times a thing was dominated to 0
        thisWasDominatedXTimes[h] = 0;
    }
    for (int i = 0; i < newAdults.size()-1; i++){ //goes through all the adults seeing which dominates
        for(int j = i+1; j < newAdults.size(); j++){ //j starts at the index past i so things are not checked against themselves and things are not rechecked
            //since a pair of indices should only be compared to one another once, it tracks the domination for both i and j
            if (dominationCheckTest(newAdults[i], newAdults[j])){ //i dominates
                arrayOfDominations[i].push_back(j);
                thisWasDominatedXTimes[j]++;
                cout << i << " dominated " << j << endl;
            }
            else if (dominationCheckTest(newAdults[j], newAdults[i])){//j dominates
                arrayOfDominations[j].push_back(i);
                thisWasDominatedXTimes[i]++;
                cout << j << " dominated " << i << endl;
            }
        }
    }
    for (int i = 0; i < genSize; i++){ //switch thisWasDominatedXTimes to the net number of times this was dominated by other things (negative numbers mean it did more dominating than getting dominated)
        thisWasDominatedXTimes[i] -= arrayOfDominations[i].size();
    }
    int leastDominated = 1000; //starts at a random value that is far larger than anything that should occur in my code - as the name implies, this represents the smallest number of times an adult was dominated
    bool all1000 = false; //when an adult has been given its permanent rank, the number of times it was dominated will be set to 1000 to ensure everything gets a rank
    int k = 0; //a number k that represents the rank that will be assigned to an individual (k is a bad name for it, that was a holdover from when the below was a for loop not a while loop)
    while (k < genSize && all1000 == false){ //loops until either everything has been given a rank
        //resets leastDominated to a bad value and all1000 to true so the loop will repeat in the same manner as the previous run
        leastDominated = 1000; 
        all1000 = true;
        for(int l = 0; l < genSize; l++){ //goes through every member of a generation, seeing how many times it was dominated
            if (k == 0){
                cout << l << ": " << thisWasDominatedXTimes[l] << endl;
            }
            if (thisWasDominatedXTimes[l] < leastDominated){ //if this was dominated fewer times than leastDominated, it is given the rank k+1 and the number of times it was dominated is the new leastDominated
                leastDominated = thisWasDominatedXTimes[l];
                all1000 = false;
                newAdults[l].rank = k+1;
            }
            else if (thisWasDominatedXTimes[l] == leastDominated && leastDominated != 1000){ //if this number was dominated the same number of times as the next least dominated, they should have the same rank (the !=1000 keeps it from reranking ones that already have a rank)
                all1000 = false;
                newAdults[l].rank = k+1;
            }
        }
        for (int anotherLoop = 0; anotherLoop < genSize; anotherLoop++){ //once the true value for leastDominated was selected, it sets any numbers who were dominated leastDomination times to saying they were dominated 1000 times so they are not reranked next 
            if (thisWasDominatedXTimes[anotherLoop] == leastDominated && leastDominated != 1000){
                thisWasDominatedXTimes[anotherLoop] = 1000;
            }
        }
        k++; //increments k so not stuck in an infinite loop
    }
    
    //prints everything so we can set it 
    for (int m = 0; m < newAdults.size(); m++){
        cout << m << ":" << newAdults[m].getRank() << " ~ ";
        if (newAdults[m].getRank() == INT_MAX){
            cout << "\nERROR - the function didn't work" << std::endl;
        }
    }
    cout << std::endl; 
}

void sortaGiveDistance(std::vector<Adult> & pool){
    int poolSize = pool.size();
    
    //starting rankSort to make sure nans are at the end of the array.
    std::sort(pool.begin(), pool.begin() + poolSize, rankSort);

    for (int i = 0; i < poolSize; i++ ){
        //reset each individual's distance
        pool[i].distance = 0.0;
    }
    
    //Sort by the first objective function, posDiff
    std::sort(pool.begin(), pool.begin() + poolSize, LowerPosDiff);
    //Set the boundaries
    pool[0].distance = 1e+12;
    pool[poolSize - 1].distance = 1e+12;


    //For each individual besides the upper and lower bounds, make their distance equal to
    //the current distance + the absolute normalized difference in the function values of two adjacent individuals.
    double normalPosDiffLeft;
    double normalPosDiffRight;
    for(int i = 1; i < poolSize - 1; i++){
        //Divide left and right individuals by the worst individual to normalize
        normalPosDiffLeft = pool[i+1].posDiff/pool[poolSize - 1].posDiff;
        normalPosDiffRight = pool[i-1].posDiff/pool[poolSize - 1].posDiff;
        //distance = distance + abs((i+1) - (i-1))
        pool[i].distance = pool[i].distance + abs((normalPosDiffLeft - normalPosDiffRight));// /(pool[poolSize - 1].posDiff - pool[0].posDiff));
    }

    //Repeat above process for speedDiff    
    std::sort(pool.begin(), pool.begin() + poolSize, LowerSpeedDiff);
    //Set the boundaries
    pool[0].distance = 1e+12;
    pool[poolSize - 1].distance = 1e+12;

    
    //For each individual besides the upper and lower bounds, make their distance equal to
    //the current distance + the absolute normalized difference in the function values of two adjacent individuals.
    double normalSpeedDiffLeft;
    double normalSpeedDiffRight;
    for(int i = 1; i < poolSize - 1; i++){
        //Divide left and right individuals by the worst individual to normalize
        normalSpeedDiffLeft = pool[i+1].speedDiff/pool[poolSize - 1].speedDiff;
        normalSpeedDiffRight = pool[i-1].speedDiff/pool[poolSize - 1].speedDiff;
        //distance = distance + abs((i+1) - (i-1))
        pool[i].distance = pool[i].distance + abs((normalSpeedDiffLeft - normalSpeedDiffRight));// /(pool[poolSize - 1].speedDiff - pool[0].speedDiff));
    }

    for (int i = 0; i < poolSize; i++){
        cout << i << ":" << pool[i].distance << " ~ ";
    }
    cout << endl;
}

bool createMasks(std::mt19937_64& rng, bool printMask){
    bool wholeRandGood = true, averageGood = true, bundleVarsGood = true, allGood = true;
    //Create a mask to determine which of the child's parameters will be inherited from which parents
    //This is done the same way it is in newGeneration in ga_crossover.cpp
    std::vector<int> mask;
    for (int j = 0; j < OPTIM_VARS; j++) //initializes mask to average in case it gets stuck - also to get the correct length preset
    {
        mask.push_back(AVG); 
    }

    crossOver_wholeRandom(mask, rng);
    for (int i = 0; i < OPTIM_VARS; i++){
        if (mask[i] != 1 && mask[i] != 2){
            wholeRandGood = false;
            allGood = false;
        }
    }
    if (!wholeRandGood){
        cout << "crossOver_wholeRandom generates the mask incorrectly" << endl;
    }

    crossOver_average(mask);
    for (int i = 0; i < OPTIM_VARS; i++){
        if (mask[i] != 3){
            averageGood = false;
            allGood = false;
        }
    }
    if (!averageGood){
        cout << "crossOver_average generates the mask incorrectly" << endl;
    }

    crossOver_bundleVars(mask, rng);
    for (int i = 0; i < OPTIM_VARS; i++){
        if (mask[i] != 1 && mask[i] != 2){
            cout << mask[i] << endl;
            bundleVarsGood = false;
            allGood = false;
        }
        if(i > GAMMA_OFFSET && i < GAMMA_OFFSET + GAMMA_ARRAY_SIZE){
            if (mask[i] != mask[i-1]){
                cout << i << ": " << mask[i] << " vs " << i-1 << ": " << mask[i-1] << endl;
                bundleVarsGood = false;
                allGood = false;
            }
        }
        else if(i > TAU_OFFSET && i < TAU_OFFSET + TAU_ARRAY_SIZE){
            if (mask[i] != mask[i-1]){
                cout << i << ": " << mask[i] << " vs " << i-1 << ": " << mask[i-1] << endl;
                bundleVarsGood = false;
                allGood = false;
            }
        }
        else if(i > COAST_OFFSET && i < COAST_OFFSET + COAST_ARRAY_SIZE){
            if (mask[i] != mask[i-1]){
                cout << i << ": " << mask[i] << " vs " << i-1 << ": " << mask[i-1] << endl;
                bundleVarsGood = false;
                allGood = false;
            }
        }
    }
    if(!bundleVarsGood){
        cout << "crossOver_bundleVars generates the mask incorrectly" << endl;
    }

    return allGood;
}

bool makeChildrenWithDifferentMethods(std::mt19937_64& rng, cudaConstants * utcConstants){

    const int expectedNumChildren = 6; //pairs of children should be generated 3 times, so the expected number of children that should be created are 6
    int childrenCreated = 0;
    //selected values for alpha (indices 0 & 1), beta (2 & 3), zeta (4 & 5), and tripTime that are within the acceptable range for rkParameters
    //these values are near the edge of the acceptable range of values to make them easier to calculate
    double parentVals[8] = {3.1414, -3.1414, 3.1414, 0.0, 1.5707, -1.5707,45000000.0, 35000000.0};
    rkParameters<double> p1(parentVals[6], parentVals[0], parentVals[2], parentVals[4]); //using unit testing constructor with only tripTime, alpha, beta, and zeta
    rkParameters<double> p2(parentVals[7], parentVals[1], parentVals[3], parentVals[5]);
    Child par1(p1);
    Child par2(p2);
    Adult parent1(par1);
    Adult parent2(par2);
    Child* children = new Child[expectedNumChildren]; 

    //sets the generation to 1 because generateChildPair takes in a generation
    // 1 was chosen because this is the first generation and the generation number should not have a major impact on the code
    int gen = 1; 

    //while there have not been any errors detected in the code this is set to true
    bool noErrors = true;

    for (int i = 0; i < 4; i++){
        cout << "i = " << i << endl;
        if (i == 2){
            utcConstants->thruster_type = 1;
        }
        if(i %2 == 1){
            utcConstants->mutation_rate = 1.0;
        }
        else{
            utcConstants->mutation_rate = 0.0;
        }
        childrenCreated = 0;

        //Create a mask to determine which of the child's parameters will be inherited from which parents
        std::vector<int> mask;
        for (int j = 0; j < OPTIM_VARS; j++) //initializes mask to average in case it gets stuck - also to get the correct length preset
        {
            mask.push_back(AVG); 
        }

        double annealing = 0.10; //the initial annealing value from the config

        //Generate a pair of children based on the random cross over mask method from a pair of parents
        //Generate the base mask
        crossOver_wholeRandom(mask, rng);

        //Generate a pair of children based on the mask
        generateChildrenPair(parent1, parent2, children, mask, annealing, rng, childrenCreated, gen, utcConstants);

        if (!checkReasonability(children[childrenCreated -2], children[childrenCreated-1], mask, parentVals, 1)){
            noErrors = false;
        }

        //Generate a pair of children from a pair of parents based on the average mask method
        //Generate the base mask
        crossOver_average(mask);

        //Generate a pair of children based on the mask
        generateChildrenPair(parent1, parent2, children, mask, annealing, rng, childrenCreated, gen, utcConstants);

       if (!checkReasonability(children[childrenCreated -2], children[childrenCreated-1], mask, parentVals, 2)){
            noErrors = false;
        }
        //Generate a pair of children from a pair of parents based on the bundled random mask method
        //Generate the base mask
        crossOver_bundleVars(mask, rng);

        //Generate a pair of children based on the mask
        generateChildrenPair(parent1, parent2, children, mask, annealing, rng, childrenCreated, gen, utcConstants);
        
        if (!checkReasonability(children[childrenCreated -2], children[childrenCreated-1], mask, parentVals, 3)){
            noErrors = false;
        }

        if((childrenCreated != expectedNumChildren && i < 2) || !noErrors){
            cout << "Could not even generate children correctly with no thrust" << endl;
            return false;
        }

        delete[] children;
    }

    if(childrenCreated == expectedNumChildren && noErrors){
        return true;
    }
    else{
        cout << "Could not generate children correctly when there was thrust" << endl;
        return false;
    }

}

//1 stands for crossOver_wholeRandom, 2 stands for crossOver_average, 3 stands for crossOver_bundleVars
bool checkReasonability(const Child& c1, const Child& c2, std::vector<int> & mask, double parentsValues[], int whichMethod){
    bool noErrors = true;
    bool skipPrint = false;
    int parValIndex = 0;
    for (int i = ALPHA_OFFSET; i <= TRIPTIME_OFFSET; i++){
        if (mask[i] == PARTNER1 && whichMethod != 2){ //the mask was lipped, so compare the second child to it first, and the other child should have the other parent's parameters
            if (c2.startParams.tripTime != parentsValues[parValIndex] || c1.startParams.tripTime != parentsValues[parValIndex+1]){
                cout << parentsValues[parValIndex] << " should equal c2 " << c2.startParams.tripTime << endl;
                cout << "Error with ";
                noErrors = false;
            }
            else{
                cout << "Expected values for ";
            }
        }
        else if (mask[i] == PARTNER2 && whichMethod != 2){
            if (c2.startParams.tripTime != parentsValues[parValIndex+1] || c1.startParams.tripTime != parentsValues[parValIndex]){
                cout << parentsValues[parValIndex+1] << " should equal c2 " << c2.startParams.tripTime << endl;
                cout << "Error with ";
                noErrors = false;
            }
            else{
                cout << "Expected values for ";
            }
        }
        else if (mask[i] == AVG && whichMethod != 1){
            if (c2.startParams.tripTime != (parentsValues[parValIndex]+parentsValues[parValIndex+1])/2 || c1.startParams.tripTime != (parentsValues[parValIndex]+parentsValues[parValIndex+1])/2){
                cout << "Error with ";
                noErrors = false;
            }
            else{
                cout << "Expected values for ";
            }
        }
        else{
            cout << "Error with mask for ";
            skipPrint = true;
        }
        if (!skipPrint){
            switch (i)
            {
            case ALPHA_OFFSET:
                cout << "alpha for ";
                break;
            case BETA_OFFSET:
                cout << "beta for ";
                break;
            case ZETA_OFFSET:
                cout << "zeta for ";
                break;
            case TRIPTIME_OFFSET:
                cout << "trip time for ";
                break;
            default:
                noErrors = false;
                break;
            }
        }
        switch (whichMethod){
        case 1:
            cout << "crossOver_wholeRandom" << endl;
            break;
        case 2:
            cout << "crossOver_average" << endl;
            break;
        case 3:
            cout << "crossOver_bundleVars" << endl;
            break;
        default:
            cout << "Error with many things including this function..." << endl;
            noErrors = false;
            break;
        }
        if (!noErrors){
            return false;
        }
        parValIndex += 2; //must increase this at the end of each time through the loop so it's comparing to the correct numbers
    }
    return true;
}

/*
bool firstFullGen(){
    Child* genZero = new Child[genSize];
    elements<double> elems(); 
    coefficients<double> coeffs();
    genZero[0] = Child(rkParameters<double>(40000000.0, elems, coeffs)); //about 1.27 years - going to call rank 3, distance 3000
    genZero[1] = Child(rkParameters<double>(35000000.0, elems, coeffs)); //about 1.12 years - going to call rank 2, distance 2700
    genZero[2] = Child(rkParameters<double>(41000000.0, elems, coeffs)); //about 1.30 years - going to call rank 3, distance 3200
    genZero[3] = Child(rkParameters<double>(32000000.0, elems, coeffs)); //about 1.01 years - going to call rank 1, distance 1000
    genZero[4] = Child(rkParameters<double>(47000000.0, elems, coeffs)); //about 1.49 years - going to call rank 4, distance 10000
    genZero[5] = Child(rkParameters<double>(43000000.0, elems, coeffs)); //about 1.36 years - going to call rank 3, distance 5000
    genZero[6] = Child(rkParameters<double>(37000000.0, elems, coeffs)); //about 1.17 years - going to call rank 2, distance 2400
    genZero[7] = Child(rkParameters<double>(38000000.0, elems, coeffs)); //about 1.20 years - going to call rank 2, distance 2300
    genZero[8] = Child(rkParameters<double>(44000000.0, elems, coeffs)); //about 1.40 years - going to call rank 4, distance 4000
    genZero[9] = Child(rkParameters<double>(45000000.0, elems, coeffs)); //about 1.43 years - going to call rank 4, distance 4500
    for (int i = 0; i < genSize; i++){
        genZero[i] = Child(paramsForIndividuals[i]);
    }
    std::vector<Adult> parents;
    UTconvertToAdult(parents, genZero);
    int ranks[genSize] = {3,2,3,1,4,3,2,2,4,4};
    int distances[genSize] = {3000,2700,3200,1000,10000,5000,2400,2300,4000,4500};
    for (int i = 0; i < genSize; i++){
        parents[i].rank = ranks[i];
        parents[i].distance = distances[i];
    }
    std::sort(parents.begin(), parents.end(), rankDistanceSort);
    std::vector<Adult> youngGen;


    delete[] genZero;
    
}
*/

bool checkUTMutateMask(){
    //set up the mask to be mutated
    bool * mutateMask= new bool[OPTIM_VARS];
    //rng based on time seed
    std::mt19937_64 rng(time(NULL));
    //mutation rate from config
    double mutation_rate = 0.75;
    //call the UT fuction
    UTmutateMask(rng, mutateMask, mutation_rate);

    //print the results
    cout << "Mask: " << endl;
    for(int i = 0; i < OPTIM_VARS; i++){
        cout << i << " M: " << mutateMask[i] << endl;
    }

    return true;
    
}

void UTmutateMask(std::mt19937_64 & rng, bool * mutateMask, double mutation_rate){
    
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
    //make sure this matched how many are set to true
    cout << "Final geneCount: " << geneCount << endl;
}