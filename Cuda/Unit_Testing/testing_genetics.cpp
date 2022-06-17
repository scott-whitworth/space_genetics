bool runGeneticsUnitTests(bool printThings){
// SETTING UP CUDA CONSTANTS TO BE USED BY OTHER FUNCTIONS 
    //making cudaConstants to control what is going on in genetic algorithm while still using the original functions
    cudaConstants* utcConstants = new cudaConstants(); 

    //sets these tolerances so that the domination in stolenGiveRank works properly
    //if these are not set, it gives the wrong rank to one Adult from each of the lower ranks (it promotes them to the next rank)
    utcConstants->posDominationTolerance = 1.0e-14;
    utcConstants->speedDominationTolerance = 1.0e-14;

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

    //Flag for what type of landing you want. Current modes: "soft"=2, "hard"=1
    //testing for soft because as of 2022 that is the main mission type we are exploring 
    //additionally, mission type should have little impact on the genetics algorithms
    utcConstants->missionType = 2; 

    //STUFF SO RUNGE_KUTTACUDA.CU CAN RUN
    utcConstants->thread_block_size =32;
    utcConstants->orbitalPeriod = 3.772645011085093e+07; // orbital period time of 1999RQ36 Bennu (s) 
    utcConstants->GuessMaxPossibleSteps = 1000000;

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
    if (firstParentsTest(utcConstants, printThings)){
        cout << "PASSED: Children can be converted to adults that can be sorted" << endl;
    }
    else{
        cout << "FAILED: Children cannot be converted to adults that can be sorted" << endl;
        allWorking = false;
    }
    if (createMasks(rng, printThings)){
        cout << "PASSED: All three masks were successfully generated as expected" << endl;
    }
    else{
        cout << "FAILED: Not all the masks were successfully generated" << endl;
        allWorking = false;
    }
    if (makeChildrenWithDifferentMethods(rng, utcConstants, printThings)){
        cout << "PASSED: Successfully made children using the three different methods" << endl;
    }
    else{
        cout << "FAILED: Could not successfully make children using the three different methods" << endl;
        allWorking = false;
    }
    if (firstFullGen(rng, utcConstants, printThings)){
        cout << "PASSED: Successfully made the first generation of adults from another set of adults" << endl;
    }
    else{
        cout << "FAILED: Could not successfully make the first generation of adults from another set of adults" << endl;
        allWorking = false;
    }


    delete[] utcConstants;
    delete launchCon;
    return allWorking;
}

//returns true if the first generation is generated
bool firstParentsTest(const cudaConstants * utcConstants, bool printThings){
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

    for (int i = 0; i < utcConstants->num_individuals; i++){
        parents.push_back(Adult(theChildren[i]));
    }

    //Ranks based on my hand calculation of ranks based on the version of giveRank from 6/14/22 in the eveing 
    //Dominations are based on the dominations test
    //Below is a list of the parameters of the above children and which rank they should belong to based on their position and speed differences
    //The work is shown in a PDF on Teams in the 2022 folder
    //Rank 1: (150, 20) distance- 1.176471; (20, 70) distance- 1.051546; (110, 40) distance - 0.247119
    //Rank 2:  (180, 30) distance - 0.226501; (30, 90) distance - 0.108854; (20, 120) distance - 0.060340 
    //Rank 3: (120, 90) distance - 1.111583; (50,120) distance - 0.117647
    //Rank 4: (220,970) distance - 1.470588; (340,90) distance - 1.030928

    stolenGiveRank(parents, utcConstants);
    //wrongWayToRank(parents);
    std::sort(parents.begin(), parents.end(), rankSort);
    stolenGiveDistance(parents, utcConstants);
    std::sort(parents.begin(), parents.end(), rankDistanceSort);

    //The final result should be... (see the PDF for the work done to calculate this)
    //Rank 1: (150, 20); (20, 70); (110, 40)
    //Rank 2:  (180, 30); (30, 90); (20, 120)  
    //Rank 3: (50,120); (120, 90)
    //Rank 4: (220,970); (340,90)
    
    if (printThings){
        for (int i = 0; i < parents.size(); i++){
            cout << "(" << parents[i].posDiff << "," << parents[i].speedDiff <<"): " <<parents[i].unitTestingRankDistanceStatusPrint() << " / ";
        }
        cout << endl;
    }

    delete[] theChildren;

    if (parents.size() == utcConstants->num_individuals && checkParentsTest(parents)){
        return true;
    }
    else{
        cout << "Could not turn children into parents to create a parent generation" << endl;
        return false;
    }

}

bool checkParentsTest(std::vector<Adult>& theResults){
    std::vector<Adult> theAnswers;
    theAnswers.push_back(Adult(Child(150, 20)));
    theAnswers.push_back(Adult(Child(20, 70)));
    theAnswers.push_back(Adult(Child(110, 40)));
    theAnswers.push_back(Adult(Child(180, 30)));
    theAnswers.push_back(Adult(Child(30, 90)));
    theAnswers.push_back(Adult(Child(20, 120)));
    theAnswers.push_back(Adult(Child(50,120)));
    theAnswers.push_back(Adult(Child(120, 90)));    
    theAnswers.push_back(Adult(Child(220,970)));
    theAnswers.push_back(Adult(Child(340,90)));


    if (theAnswers.size() != theResults.size()){
        return false;
    }

    for (int i = 0; i < theResults.size(); i++){
        if (theResults[i].posDiff != theAnswers[i].posDiff || theResults[i].speedDiff != theAnswers[i].speedDiff){
            cout << "The Adult here is (" << theResults[i].posDiff << "," << theResults[i].speedDiff << "), but it should be (" << theAnswers[i].posDiff << "," << theAnswers[i].speedDiff << ")" << endl;
            return false;
        }
    }
    
    return true;
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
        if (printMask){
            cout << mask[i] << " ";
        }
        if (mask[i] != 1 && mask[i] != 2){
            wholeRandGood = false;
            allGood = false;
        }
    }
    if (!wholeRandGood){
        cout << "\ncrossOver_wholeRandom generates the mask incorrectly" << endl;
    }
    else if (printMask){
        cout << endl;
    }

    crossOver_average(mask);
    for (int i = 0; i < OPTIM_VARS; i++){
        if (printMask){
            cout << mask[i] << " ";
        }
        if (mask[i] != 3){
            averageGood = false;
            allGood = false;
        }
    }
    if (!averageGood){
        cout << "crossOver_average generates the mask incorrectly" << endl;
    }
    else if (printMask){
        cout << endl;
    }

    crossOver_bundleVars(mask, rng);
    for (int i = 0; i < OPTIM_VARS; i++){
        if (printMask){
            cout << mask[i] << " ";
        }
        if (mask[i] != 1 && mask[i] != 2){
            bundleVarsGood = false;
            allGood = false;
        }
        if(i > GAMMA_OFFSET && i < GAMMA_OFFSET + GAMMA_ARRAY_SIZE){
            if (mask[i] != mask[i-1]){
                cout << "\n" << i << ": " << mask[i] << " vs " << i-1 << ": " << mask[i-1] << endl;
                bundleVarsGood = false;
                allGood = false;
            }
        }
        else if(i > TAU_OFFSET && i < TAU_OFFSET + TAU_ARRAY_SIZE){
            if (mask[i] != mask[i-1]){
                cout << "\n" << i << ": " << mask[i] << " vs " << i-1 << ": " << mask[i-1] << endl;
                bundleVarsGood = false;
                allGood = false;
            }
        }
        else if(i > COAST_OFFSET && i < COAST_OFFSET + COAST_ARRAY_SIZE){
            if (mask[i] != mask[i-1]){
                cout << "\n" << i << ": " << mask[i] << " vs " << i-1 << ": " << mask[i-1] << endl;
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

bool makeChildrenWithDifferentMethods(std::mt19937_64& rng, cudaConstants * utcConstants, bool printThings){
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

    //while there have not been any errors detected in the codel this is set to true
    bool noErrors = true;

    //these individuals should not be mutated or else we cannot verify that generating individuals is working as expected
    utcConstants->mutation_rate = 0.0; 

    for (int i = 0; i < 2; i++){
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

        //the "magic" 1 represents the first test, crossover_wholeRandom
        if (!checkReasonability(children[childrenCreated -2], children[childrenCreated-1], mask, parentVals, 1, utcConstants, printThings)){
            noErrors = false;
        }

        //Generate a pair of children from a pair of parents based on the average mask method
        //Generate the base mask
        crossOver_average(mask);

        //Generate a pair of children based on the mask
        generateChildrenPair(parent1, parent2, children, mask, annealing, rng, childrenCreated, gen, utcConstants);

        //the "magic" 2 represents the second test, crossover_average
       if (!checkReasonability(children[childrenCreated -2], children[childrenCreated-1], mask, parentVals, 2, utcConstants, printThings)){
            noErrors = false;
        }
        //Generate a pair of children from a pair of parents based on the bundled random mask method
        //Generate the base mask
        crossOver_bundleVars(mask, rng);

        //Generate a pair of children based on the mask
        generateChildrenPair(parent1, parent2, children, mask, annealing, rng, childrenCreated, gen, utcConstants);
        
        //the "magic" 3 represents the third test, crossover_bundleVars
        if (!checkReasonability(children[childrenCreated -2], children[childrenCreated-1], mask, parentVals, 3, utcConstants, printThings)){
            noErrors = false;
        }

        //if there were errors or not enough children were created on the first run with no thruster, won't bother testing with thruster
        //the thruster should not make any difference in the parameters being checked by checkReasonability
        if((childrenCreated != expectedNumChildren && i == 0) || !noErrors){
            cout << "Could not even generate children correctly with no thruster" << endl;
            delete[] children;
            return false;
        }
        
        //once it has gone through once with no thruster, it turns the thuster on to make sure it works with the thruster too
        utcConstants->thruster_type = 1; 
    }

    //if it executes correctly, deallocates dynamic memory and returns true
    if(childrenCreated == expectedNumChildren && noErrors){
        delete[] children;
        return true;
    }
    //if there were errors, it deallocates dynamic memory and returns false
    else{
        cout << "Could not generate children correctly when there was thrust" << endl;
        delete[] children;
        return false;
    }

}

//1 stands for crossOver_wholeRandom, 2 stands for crossOver_average, 3 stands for crossOver_bundleVars
bool checkReasonability(const Child& c1, const Child& c2, std::vector<int> & mask, double parentsValues[], int whichMethod, cudaConstants * utcConstants, bool printThings){
    bool noErrors = true; //initially there are no error detected yet
    bool notFirstSwitch = false; //this is just here so the output turns out nicely
    bool skipPrint = !printThings; //skipping printing is the opposite of printing
    int parValIndex = 0;
    for (int i = ALPHA_OFFSET; i <= TRIPTIME_OFFSET; i++){
        //the mask was flipped, so compare the second child to it first, and the other child should have the other parent's parameters
        //because flipping a mask sets PARTNER1 values to PARTNER2 and vice versa
        if (mask[i] == PARTNER1 && whichMethod != 2){ 
            //if the mask at this index corresponds to PARTNER1, the second child generated should have parent 1's value for this parameter
            //and the first child should have the values the other parent (PARTNER2) has for this parameter
            //if either of these are not the case, there is an issue
            if (getParamStuff(i,c2) != parentsValues[parValIndex] || getParamStuff(i,c1) != parentsValues[parValIndex+1]){ 
                cout << parentsValues[parValIndex] << " should equal c2 " << getParamStuff(i,c2) << endl;
                cout << "Error with ";
                skipPrint = false; //whether or not valid messages are being printed, it needs to print error message
                noErrors = false;
            }
            else if (printThings){
                cout << "Expected values for ";   
            }
            else{
                skipPrint = true;
            }
        }
        else if (mask[i] == PARTNER2 && whichMethod != 2){
            //if the mask at this index corresponds to PARTNER2, the second child generated should have parent 2's value for this parameter
            //and the first child should have the values the other parent (PARTNER1) has for this parameter
            //if either of these are not the case, there is an issue
            if (getParamStuff(i,c2) != parentsValues[parValIndex+1] || getParamStuff(i,c1) != parentsValues[parValIndex]){
                cout << parentsValues[parValIndex+1] << " should equal c2 " << getParamStuff(i,c2) << endl;
                cout << "Error with ";
                skipPrint = false; //whether or not valid messages are being printed, it needs to print error message
                noErrors = false;
            }
            else if (printThings){
                cout << "Expected values for ";   
            }
            else{
                skipPrint = true;
            }
        }
        //the only time it would not be an error for a child to have AVG as a value in its mask is when the method being used is crossOver_average (whichMethod == 2)
        else if (mask[i] == AVG && whichMethod == 2){
            //if the value is supposed to be an average, ensures it actually is
            if (getParamStuff(i,c2) != (parentsValues[parValIndex]+parentsValues[parValIndex+1])/2 || getParamStuff(i,c1) != (parentsValues[parValIndex]+parentsValues[parValIndex+1])/2){
                cout << "Error with ";
                skipPrint = false; //whether or not valid messages are being printed, it needs to print error message
                noErrors = false;
            }
            else if (printThings){
                cout << "Expected values for ";   
            }
            else{
                skipPrint = true;
            }
        }
        //if the values in the mask does not correspond to PARTNER1,PARTNER2, or AVG then there is an issue with the mask
        else{
            cout << "Error with mask for ";
            noErrors = false;
            skipPrint = false; //whether or not valid messages are being printed, it needs to print error message
            //it is not an issue with the values for alpha, beta, zeta, or tripTime necessarily, so it will not print the cout about the type
            notFirstSwitch = true; 
        }
        if (!skipPrint){ //if a message should be printed, goes through these statements
            if (!notFirstSwitch){
                switch (i)
                {
                case ALPHA_OFFSET:
                    cout << "alpha corresponded with the actual values for ";
                    break;
                case BETA_OFFSET:
                    cout << "beta corresponded with the actual values for ";
                    break;
                case ZETA_OFFSET:
                    cout << "zeta corresponded with the actual values for ";
                    break;
                case TRIPTIME_OFFSET:
                    cout << "trip time corresponded with the actual values for ";
                    break;
                default:
                    cout << "issue with verification function... ";
                    noErrors = false; //something weird is going on with i which is also an error
                    break;
                }
            }
            switch (whichMethod){
            case 1:
                cout << "crossOver_wholeRandom";
                break;
            case 2:
                cout << "crossOver_average";
                break;
            case 3:
                cout << "crossOver_bundleVars";
                break;
            default:
                cout << "incorrectly specified cross over method";
                noErrors = false; //if the method was not correctly specified, this is an error
                break;
            }
            if (utcConstants->thruster_type == thruster<double>::NO_THRUST){
                cout << " with no thruster" <<endl;
            }
            else {
                cout << " with a thruster" << endl;
            }
        }
        if (!noErrors){ //if there were errors, does not bother to loop again because likely the other parameters will be broken 
            return false;
        }
        parValIndex += 2; //must increase this at the end of each time through the loop so it's comparing to the correct numbers
    }
    return true;
}

//just makes checkReasonability shorter - takes in an offset and accesses a child's parameter corresponding to this offset
double getParamStuff(const int correspondingOffset, const Child& aChild){
    if (correspondingOffset == ALPHA_OFFSET){ 
        return aChild.startParams.alpha;
    }
    else if (correspondingOffset == BETA_OFFSET){
        return aChild.startParams.beta;
    }
    else if (correspondingOffset == ZETA_OFFSET){
        return aChild.startParams.zeta;
    }
    else if (correspondingOffset == TRIPTIME_OFFSET){
        return aChild.startParams.tripTime;
    }
    else { //if there is an error with what is passed in, returns an error message and a weird, small negative
        cout << "ERROR: undefined offset" << endl;
        return -1e-06; //there was no specific reason for this number and it can be changed
    }
}

//unit testing creating an entire generation of individuals
bool firstFullGen(std::mt19937_64& rng, cudaConstants * utcConstants, bool printThings){
    bool noErrors = true;

    Child* genZero = new Child[utcConstants->num_individuals]; 
    elements<double> elems; //makes a set of default elements
    coefficients<double> coeffs = coefficients<double>(); //makes a set of default coefficients

    rkParameters<double>* paramsForIndividuals = new rkParameters<double>[utcConstants->num_individuals];
    
    //creates parameters for making child with unique tripTimes, but default elements and default coefficients 
    paramsForIndividuals[0] = rkParameters<double>(40000000.0, elems, coeffs); //about 1.27 years  
    paramsForIndividuals[1] = rkParameters<double>(35000000.0, elems, coeffs); //about 1.12 years 
    paramsForIndividuals[2] = rkParameters<double>(41000000.0, elems, coeffs); //about 1.30 years
    paramsForIndividuals[3] = rkParameters<double>(32000000.0, elems, coeffs); //about 1.01 years 
    paramsForIndividuals[4] = rkParameters<double>(47000000.0, elems, coeffs); //about 1.49 years 
    paramsForIndividuals[5] = rkParameters<double>(43000000.0, elems, coeffs); //about 1.36 years 
    paramsForIndividuals[6] = rkParameters<double>(37400000.0, elems, coeffs); //about 1.17 years 
    paramsForIndividuals[7] = rkParameters<double>(38000000.0, elems, coeffs); //about 1.20 years 
    paramsForIndividuals[8] = rkParameters<double>(44000000.0, elems, coeffs); //about 1.40 years 
    paramsForIndividuals[9] = rkParameters<double>(45000000.0, elems, coeffs); //about 1.43 years 
    
    //turns the rkParameters into children to make up genZerp
    for (int i = 0; i < utcConstants->num_individuals; i++){
        genZero[i] = Child(paramsForIndividuals[i], utcConstants);
    }

    std::vector<Adult> parents;
    firstGeneration(genZero, parents, utcConstants); //genZero and is turned into parents using firstGeneration which has been unit tested previously

    //the parents are sprted so they can be selected to be parents for a new generation
    stolenGiveRank(parents, utcConstants); 
    std::sort(parents.begin(), parents.end(), rankSort); 
    stolenGiveDistance(parents, utcConstants); 
    std::sort(parents.begin(), parents.end(), rankDistanceSort);

    std::vector<Adult> youngGen;

    //creates young a young generation with an annealing rate of 0.0 and the generation is set to 1 (the random 0.0 and 1)
    //neither of these numbers should really affect the newGeneration being created, because mutation rate and annealing are both set to 0 
    newGeneration(parents, youngGen, 0.0, 1, rng, utcConstants);

    //verifies 
    noErrors = verifyFullGen(youngGen, parents, utcConstants, printThings);
    if (youngGen.size() != parents.size() || parents.size() != utcConstants->num_individuals){
        noErrors = false;
    }

    //TODO: create a function that verifies this is working correctly  

    //deallocates dynamic memory
    delete[] genZero;
    delete[] paramsForIndividuals;

    makeManyChildren(rng, youngGen, parents, utcConstants, printThings);

    //returns noErrors
    return noErrors;
}

// A function that is used to verify that firstFullGen is working correctly
bool verifyFullGen(std::vector<Adult>& youngGen, std::vector<Adult>& possParents, const cudaConstants * utcConstants, bool printThings){
    bool noErrors = true;

    std::vector<int> parentIndexSets; //vector that holds the indices of sets of parents

    int childrenPerSet = 6;
    //int setsOfChildren = (utcConstants->num_individuals)/childrenPerSet;
    int setNum = 0;

    int firstAvgInd = 2;
    int indexVecSize = 0;

    //it checks which parents are being used and prints things if printThings is set to true
    for (int i = 0; i < utcConstants->num_individuals; i++){
        if (printThings){
            if (i < utcConstants->survivor_count){
                cout << "*|";
            }
            cout << "parents[" << i << "]'s tripTime: ";
            cout << possParents[i].startParams.tripTime;
            if (i < utcConstants->survivor_count){
                cout << "|*";
            }
            else {
                cout << "\t";
            }
            cout << "\t\t";

            cout << "youngGen[" << i << "]'s tripTime: ";
            cout << youngGen[i].startParams.tripTime << endl;
        }
        if (i == firstAvgInd + childrenPerSet*setNum){
            indexVecSize += 2;
            if (youngGen[i].startParams.tripTime == youngGen[6*setNum].startParams.tripTime){
                cout << "The parent did not have a partner and probably created clones" << endl;
                noErrors = false;
                for (int j = 0; j < utcConstants->survivor_count; j++){
                    if (youngGen[i].startParams.tripTime == possParents[j].startParams.tripTime){
                        cout << "Set " << setNum << " created by parent[" << j << "] and parent[" << j << "] " << endl;
                        parentIndexSets.push_back(j);
                        parentIndexSets.push_back(j);
                    }
                }
            }
            else{
                for (int j = 0; j < utcConstants->survivor_count; j++){
                    for (int k = j+1; k < utcConstants->survivor_count; k++){
                        if (youngGen[i].startParams.tripTime == (possParents[j].startParams.tripTime + possParents[k].startParams.tripTime)/2){
                            cout << "Set " << setNum << " created by parent[" << j << "] and parent[" << k << "] " << endl;
                            parentIndexSets.push_back(j);
                            parentIndexSets.push_back(k);
                            break;
                        }
                    }
                }
            }
            if(indexVecSize != parentIndexSets.size()){
                noErrors = false;
                cout << indexVecSize << " vs " << parentIndexSets.size() << endl;
                cout << "Unknown parents for set #" << setNum << " of children" << endl;

            }
        }
        if (i%childrenPerSet == 0){
            setNum++;
        }
    }
    return noErrors;
}

bool makeManyChildren(std::mt19937_64& rng, std::vector<Adult>& youngGen, std::vector<Adult>& possParents, cudaConstants * utcConstants, bool printThings){
    bool noErrors = true; 

    utcConstants->num_individuals = 65;
    utcConstants->survivor_count = 7;

    newGeneration(possParents, youngGen, 0.0, 1, rng, utcConstants);
    //verifies 
    noErrors = verifyFullGen(youngGen, possParents, utcConstants, false);

    if (!noErrors){
        cout << "Problem with 65 individuals generated from 7 parents" << endl;
    }

    utcConstants->survivor_count = 5;
    newGeneration(possParents, youngGen, 0.0, 1, rng, utcConstants);
    //verifies 
    noErrors = verifyFullGen(youngGen, possParents, utcConstants, false);

    if (!noErrors){
        cout << "Problem with 65 individuals generated from 5 parents" << endl;
    }

    return noErrors;
}

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

//took this directly from optimization.cu on 6/16/22 at 11:09am - will become out of date if changes made to version in optimization.cu
void stolenGiveRank(std::vector<Adult> & allAdults, const cudaConstants* cConstants) {
    //non-denominated sorting method
    //https://www.iitk.ac.in/kangal/Deb_NSGA-II.pdf

    //Used to store the current front of adults. first filled with the first front adults(best out of all population)
    // filled with index of adults in allAdults
    std::vector<int> front;

    // probably a 2D vector, or an array of vectors
    //This 2D vector will store which other adults each adult has dominated
    //1st dimension will be a spot for each adult in allAdults
    //2nd dimension will store the indexes of adults in allAdults that the adult in the 1st dimension has dominated
    std::vector<std::vector<int>> domination; 
    domination.resize(allAdults.size());

    //This vector will keep track of how many times each adult in oldAdults has been dominated by another adult
    //Each index in this vector will correspond to the same index within allAdults
    std::vector<int> dominatedByCount;
    
    //fill the vector with 0s to make sure the count is accurate
    dominatedByCount.resize(allAdults.size(), 0);

    //loop through each individual within the allAdults vector
    for (int i = 0; i < allAdults.size(); i++){

        //For each individual within allAdults, compare them to each other adult
        for(int j = 0; j < allAdults.size(); j++){

            //check the status of both i and j and see if i is automatically dominated
            if(allAdults[i].errorStatus != VALID && allAdults[j].errorStatus == VALID){
                dominatedByCount[i]++;

            }//check the status of both i and j and see if j is automatically dominated
            else if(allAdults[j].errorStatus != VALID && allAdults[i].errorStatus == VALID){
                domination[i].push_back(j);
                
            }//if either both are valid or both are not valid, it will rank them normally
            //Check to see if i dominates j
            else if (dominationCheck(allAdults[i], allAdults[j], cConstants)){
                //Put the jth index in the set of individuals dominated by i
                //std::cout << "\n" << i << "th (i) Adult dominates " << j << "th (j) Adult!\n";
                domination[i].push_back(j);
            }
            //Check to see if j dominates i
            else if ( dominationCheck( allAdults[j], allAdults[i], cConstants) ){
                //std::cout << "\n" << j << "th (j) Adult dominates " << i << "th (i) Adult!\n";
                //Add one to i's dominated by count
                dominatedByCount[i]++; 
            }
        }
        
        //if i was never dominated, add it's index to the best front, front1. Making its ranking = 1.
        if (dominatedByCount[i] == 0){
            allAdults[i].rank = 1;
            front.push_back(i);

            //std::cout << "\n\nAdult #" << i << " ranked " << 1;
        }
    }

    //Used to assign rank number
    int rankNum = 1;
    //vector to store individuals' indexes in next front
    std::vector<int> newFront;

    //go until all individuals have been put in better ranks and none are left to give a ranking
    while(!front.empty()) {
        //empty the new front to put new individuals in
        std::vector<int>().swap(newFront);

        //loop through all individuals in old front
        //These individuals already have their rank set
        for(int k = 0; k < front.size(); k++){

            //loop through all the individuals that the individual in the old front dominated
            for(int l = 0; l < domination[front[k]].size(); l++){

                //subtract 1 from the dominated individuals' dominatedCount.
                //if an individual was dominated only once for example, it would be on the second front of individuals.
                dominatedByCount[domination[front[k]][l]]--;
                
                //if the dominated count is at 0, add the individual to the next front and make its rank equal to the next front number.
                if (dominatedByCount[domination[front[k]][l]] == 0){
                    //Assign a rank to the new most dominating adult left
                    allAdults[domination[front[k]][l]].rank = rankNum + 1;

                    //std::cout << "\n\nAdult #" << domination[front[k]][l] << " ranked " << rankNum+1;

                    //Add the index of the this adult to newFront
                    newFront.push_back(domination[front[k]][l]);                        
                }
            }
        }
        //increment the rank number
        rankNum++;
        
        //empty the current front
        std::vector<int>().swap(front);

        //Equate the current (now empty) front to the new front to transfer the indexes of the adults in newFront to front
        front = newFront;
    }
    //std::cout << "\n~~~~~~~~~~~~~~~~~~~~~~~~FINISHED RANKING~~~~~~~~~~~~~~~~~~~~~~~~\n";
}

//took this directly from optimization.cu on 6/16/22 at 11:29AM - will become out of date if changes made to version in optimization.cu
//could not directly include the one from optimization.cu because it causes errors since there is a main there and in testing_main.cpp
void stolenGiveDistance(std::vector<Adult> & allAdults, const cudaConstants* cConstants){
        //pool that holds the indexes of the valid adults
    std::vector<int> validAdults;

    //starting rankSort to make sure nans are at the end of the array.
    std::sort(allAdults.begin(), allAdults.begin() + allAdults.size(), rankSort);

    //checks if the adult is valid and then adds that index to the vector
    //the size of this vector will be used to find the distance for valid adults only
    for(int i = 0; i < allAdults.size(); i++){
        if(allAdults[i].errorStatus == VALID){
            validAdults.push_back(i);
        }
    }

    //set all to zero including the invalid adults
    for (int i = 0; i < allAdults.size(); i++ ){
        //reset each individual's distance
        allAdults[i].distance = 0.0;
    }
    
    //Sort by the first objective function, posDiff
    std::sort(allAdults.begin(), allAdults.begin() + validAdults.size(), LowerPosDiff);
    //Set the boundaries
    allAdults[0].distance += MAX_DISTANCE; //+=1
    allAdults[validAdults.size() - 1].distance += MAX_DISTANCE; //+=1


    //For each individual besides the upper and lower bounds, make their distance equal to
    //the current distance + the absolute normalized difference in the function values of two adjacent individuals.
    double normalPosDiffLeft;
    double normalPosDiffRight;
    for(int i = 1; i < validAdults.size() - 1; i++){
        //Divide left and right individuals by the worst individual to normalize
        normalPosDiffLeft = allAdults[i+1].posDiff/allAdults[validAdults.size() - 1].posDiff;
        normalPosDiffRight = allAdults[i-1].posDiff/allAdults[validAdults.size() - 1].posDiff;
        //distance = distance + abs((i+1) - (i-1))
        allAdults[i].distance = allAdults[i].distance + abs((normalPosDiffLeft - normalPosDiffRight));// /(allAdults[validAdults.size() - 1].posDiff - allAdults[0].posDiff));
    }

    //Repeat above process for speedDiff
    if(cConstants->missionType == Rendezvous){//only do this for the rendezvous mission since it has 2 objectives
        std::sort(allAdults.begin(), allAdults.begin() + validAdults.size(), LowerSpeedDiff);
        //Set the boundaries
        allAdults[0].distance += MAX_DISTANCE; //+=1
        allAdults[validAdults.size() - 1].distance += MAX_DISTANCE; //+=1

    
        //For each individual besides the upper and lower bounds, make their distance equal to
        //the current distance + the absolute normalized difference in the function values of two adjacent individuals.
        double normalSpeedDiffLeft;
        double normalSpeedDiffRight;
        for(int i = 1; i < validAdults.size() - 1; i++){
            //Divide left and right individuals by the worst individual to normalize
            normalSpeedDiffLeft = allAdults[i+1].speedDiff/allAdults[validAdults.size() - 1].speedDiff;
            normalSpeedDiffRight = allAdults[i-1].speedDiff/allAdults[validAdults.size() - 1].speedDiff;
            //distance = distance + abs((i+1) - (i-1))
            allAdults[i].distance = allAdults[i].distance + abs((normalSpeedDiffLeft - normalSpeedDiffRight));// /(allAdults[validAdults.size() - 1].speedDiff - allAdults[0].speedDiff));
        }
    }
}

/*
void wrongWayToRank(std::vector<Adult> & newAdults){
    int genSize = 10;
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
*/