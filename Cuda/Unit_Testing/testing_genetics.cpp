bool runGeneticsUnitTests(bool printThings){
// SETTING UP CUDA CONSTANTS TO BE USED BY OTHER FUNCTIONS 
    //making cudaConstants to control what is going on in genetic algorithm while still using the original functions
    cudaConstants* utcConstants = new cudaConstants(); 

    //sets these tolerances so that the domination in stolenGiveRank works properly
    //if these are not set, it gives the wrong rank to one Adult from each of the lower ranks (it promotes them to the next rank)
    //utcConstants->posDominationTolerance = 1.0e-14;
    //utcConstants->speedDominationTolerance = 1.0e-14;

    //chose pos_threshold and speed_threshold so anything past the first significant figure in the clone separation test is trivial
    //utcConstants->pos_threshold = 1.0e-11;
    //utcConstants->speed_threshold = 1.0e-11;

    // Seed used for randomization rng things, using seed 0 for consistancy / tracability 
    utcConstants->time_seed = 0; 

    //values taken directly from the genetic config file
    utcConstants->anneal_initial = 0.0; // initial value for annealing, was 0.1 in config, set to 0 because no mutations wanted at first
    utcConstants->anneal_final = 1.0e-7;   // final value for annealing, was 1.0e-7 in config, set to 0 because no mutations wanted at first

    // The percentage for probability of mutating a gene in a new individual, called iteratively to mutate more genes until the check fails
    // Starting this off at 0 because to ensure that there are no mutations initially to verify that children are generated as expected
    // Changes to 1.0 later in the code so that we can ensure mutation is working
    utcConstants->default_mutation_chance = 0; 
    utcConstants->mutation_amplitude = 1;
   
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
    //utcConstants->missionType = 2; 

    //STUFF SO RUNGE_KUTTACUDA.CU CAN RUN
    utcConstants->thread_block_size =32;
    utcConstants->orbitalPeriod = 3.772645011085093e+07; // orbital period time of 1999RQ36 Bennu (s) 
    //utcConstants->GuessMaxPossibleSteps = 1000000;

    //ALL THE STUFF SO EARTH'S INITIAL POSITION CAN BE CALCULATED
    utcConstants->triptime_max=2.0* SECONDS_IN_YEAR;
    utcConstants->triptime_min=1.0* SECONDS_IN_YEAR;
    utcConstants->timeRes=3600; // Earth Calculations Time Resolution Value
    utcConstants->v_escape = 2162.4/AU; //sqrt of DART Mission c3energy because there are no c3energy values recorded for the bennu mission
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
    //utcConstants->GuessMaxPossibleSteps=1000000; //TODO: Long term, I think we need to get rid of this, we use max_numsteps instead
    utcConstants->max_numsteps=1000;
    utcConstants->min_numsteps=400;

    //New ga_crossover stuff
    //utcConstants->distanceTolerance=1.0e-14;
    //default_mutation_factor=1.0         //Set at one to allow for the mutation scales below to be applied normally/to children at face value
    //duplicate_mutation_factor=1.0;       //Set at 10; when duplicate adults are mutated, they will undergo stronger mutations


// SETTING A RANDOM NUMBER GENERATOR (rng) TO BE USED BY FUNCTIONS
    // This rng object is used for generating all random numbers in the genetic algorithm, passed in to functions that need it
    std::mt19937_64 rng(utcConstants->time_seed);

    //VERY complicated part of the code with some possibility of errors -> just needed for the child constructor with rkParameters and cConstants as its arguments
    launchCon = new EarthInfo(utcConstants); 

// CALLING THE DIFFERENT UNIT TESTING ALGORITHMS
    bool allWorking = true;

    std::vector<Adult> aAdults;
    
    //firstParentsTest takes in the cuda constants and will verify that children can be converted into parents and sorted using rankDistanceSort
    if (firstParentsTest(utcConstants, printThings)){
        cout << "PASSED: Children can be converted to adults that can be sorted" << endl;
    }
    else{
        cout << "FAILED: Children cannot be converted to adults that can be sorted" << endl;
        allWorking = false;
    }
    //createMasks uses rng to create masks and ensure that their elements are as expected -> there are not AVGs where there should not be, etc
    if (createMasks(rng, printThings)){
        cout << "PASSED: All three masks were successfully generated as expected" << endl;
    }
    else{
        cout << "FAILED: Not all the masks were successfully generated" << endl;
        allWorking = false;
    }
    //uses rng and cuda constants to create children using different masks
    //these children are converted to Adults and then tested to ensure their values match the expected values
    if (makeChildrenWithDifferentMethods(rng, utcConstants, printThings)){
        cout << "PASSED: Successfully made children using the three different methods" << endl;
    }
    else{
        cout << "FAILED: Could not successfully make children using the three different methods" << endl;
        allWorking = false;
    }
    if (verifyChildrenFromCrossover(rng, printThings, utcConstants)){
        cout << "PASSED: Children were generated successfully using the childrenFromCrossover function" << endl;
    }
    else{
        cout << "FAILED: Children were not correctly generated using the childrenFromCrossover function" << endl;
        allWorking = false;
    }
    if (testingCopyOfNewGen(printThings, rng, utcConstants)){
        cout << "PASSED: newGeneration (without callRK) seems to be working as expected" << endl;
    }
    else{
        cout << "FAILED: newGeneration (without callRK) does not seem to be working" << endl;
    }
    
    //creates a generation of parents and then creates children from these parents
    //then these children are sent through a function that verifies their tripTime values 
    //tripTime was chosen because it was the paramerter set in creating the parents and it is easiest to see if the children have the correct values for tripTime 
    if (firstFullGen(rng, utcConstants, printThings)){
        cout << "PASSED: Successfully made the first generation of adults from another set of adults" << endl;
    }
    else{
        cout << "FAILED: Could not successfully make the first generation of adults from another set of adults" << endl;
        allWorking = false;
    }

    delete utcConstants;
    delete launchCon;
    return allWorking;
}

//returns true if the first generation is generated
bool firstParentsTest(const cudaConstants * utcConstants, bool printThings){
    std::vector<Adult> parents;
    Child* theChildren = new Child[utcConstants->num_individuals];

    //Made 10 children and set them with semi-random position difference and speed difference values
    //        that can easily be determined are either true or false
    //These numbers are semi-random multiples of 10 and were entered in no specific order at this point
    //This sort of simulates when the initial generation of oldAdults are created using random parameters 
    //(but these don't change from run to run and are easier to do compare)
    theChildren[0] = Child(1.5000, .2000); //
    theChildren[1] = Child (1.2000, .9000); //
    theChildren[2] = Child(1.8000, .3000); //
    theChildren[3] = Child (.2000, .7000); //
    theChildren[4] = Child(1.1000, .4000); //
    theChildren[5] = Child (.3000, .9000); //
    theChildren[6] = Child(.5000, 1.2000); //
    theChildren[7] = Child (2.2000, 9.7000);
    theChildren[8] = Child(.2000, 1.2000); //
    theChildren[9] = Child (3.4000, .9000); //
    for (int i = 0; i < utcConstants->num_individuals; i++){
        parents.push_back(Adult(theChildren[i]));
    }
    //Ranks based on my hand calculation of ranks based on the version of giveRank from 6/16/22 
    //Dominations are based on the dominations test
    //Below is a list of the parameters of the above children and which rank they should belong to 
    //    based on their position and speed differences
    //The work is shown in a PDF on Teams in the 2022 folder
    //Rank 1: (150, 20) distance- 1.176471; (20, 70) distance- 1.051546; (110, 40) distance - 0.247119
    //Rank 2:  (180, 30) distance - 0.226501; (30, 90) distance - 0.108854; (20, 120) distance - 0.060340 
    //Rank 3: (120, 90) distance - 1.111583; (50,120) distance - 0.117647
    //Rank 4: (220,970) distance - 1.470588; (340,90) distance - 1.030928

    giveRank(parents, utcConstants);
    giveDistance(parents, utcConstants);
    std::sort(parents.begin(), parents.end(), rankDistanceSort);
    
    //The final result should be... (see the PDF for the work done to calculate this)
    //Rank 1: (150, 20); (20, 70); (110, 40)
    //Rank 2:  (180, 30); (30, 90); (20, 120)  
    //Rank 3: (50,120); (120, 90)
    //Rank 4: (220,970); (340,90)
    
    if (printThings){
        for (int i = 0; i < parents.size(); i++){
            cout << "(" << parents[i].posDiff << "," << parents[i].speedDiff <<"): " <<parents[i].unitTestingRankDistanceStatusPrint() << " | ";
        }
        cout << endl;
    }

    delete[] theChildren;

    //Ensures the size is correct and that the parents match their expected values
    if (utcConstants->num_individuals == parents.size() && checkParentsTest(parents)){
        return true;
    }
    else{
        cout << "Could not turn children into parents to create a parent generation" << endl;
        return false;
    }

}

//A function that hold the expected values for parents test and that will compare the results
//      from firstParentsTest to the correct order of the values 
bool checkParentsTest(std::vector<Adult>& theResults){
    std::vector<Adult> theAnswers;
    theAnswers.push_back(Adult(Child(1.50, .20)));
    theAnswers.push_back(Adult(Child(.20, .70)));
    theAnswers.push_back(Adult(Child(1.10, .40)));
    theAnswers.push_back(Adult(Child(1.80, .30)));
    theAnswers.push_back(Adult(Child(.30, .90)));
    theAnswers.push_back(Adult(Child(.20, 1.20)));
    theAnswers.push_back(Adult(Child(.50,1.20)));
    theAnswers.push_back(Adult(Child(1.20, .90)));    
    theAnswers.push_back(Adult(Child(2.20,9.70)));
    theAnswers.push_back(Adult(Child(3.40,.90)));

    //if the two vectors are not the same size, we have a problem
    if (theAnswers.size() != theResults.size()){
        return false;
    }

    //loops through the results and compares its values to the expected values 
    for (int i = 0; i < theResults.size(); i++){
        //if either the posDiff or speedDiff is wrong, then the vector has been sorted wrong and it prints what the issie is
        if (theResults[i].posDiff != theAnswers[i].posDiff || theResults[i].speedDiff != theAnswers[i].speedDiff){
            cout << "The Adult here is (" << theResults[i].posDiff << "," << theResults[i].speedDiff << "), but it should be (" << theAnswers[i].posDiff << "," << theAnswers[i].speedDiff << ")" << endl;
            return false;
        }
    }
    return true;
}

//creates masks to ensure that they are generated properly
bool createMasks(std::mt19937_64& rng, bool printMask){
    bool wholeRandGood = true, averageGood = true, bundleVarsGood = true, averageRatioGood = true, allGood = true;
    //Create a mask to determine which of the child's parameters will be inherited from which parents
    //This is done the same way it is in newGeneration in ga_crossover.cpp
    std::vector<int> mask;
    for (int j = 0; j < OPTIM_VARS; j++) //initializes mask to average in case it gets stuck - also to get the correct length preset
    {
        mask.push_back(AVG); 
    }

    crossOver_wholeRandom(mask, rng);
    //everything in crossOver_wholeRandom should be either 1s or 2s
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
    //all the values in the indices in crossOver_average should contain threes
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
    //all the values in the indices in crossOver_bundleVars should contain sets of 1s or 2s
    for (int i = 0; i < OPTIM_VARS; i++){
        if (printMask){
            cout << mask[i] << " ";
        }
        //if it is neither a 1 or a 2, there is definitely an error
        if (mask[i] != 1 && mask[i] != 2){
            bundleVarsGood = false;
            allGood = false;
        }
        //gamma should all be taken from the same parent
        if(i > GAMMA_OFFSET && i < GAMMA_OFFSET + GAMMA_ARRAY_SIZE){
            if (mask[i] != mask[i-1]){
                cout << "\n" << i << ": " << mask[i] << " vs " << i-1 << ": " << mask[i-1] << endl;
                bundleVarsGood = false;
                allGood = false;
            }
        }
        //tau should all be taken from the same parent
        else if(i > TAU_OFFSET && i < TAU_OFFSET + TAU_ARRAY_SIZE){
            if (mask[i] != mask[i-1]){
                cout << "\n" << i << ": " << mask[i] << " vs " << i-1 << ": " << mask[i-1] << endl;
                bundleVarsGood = false;
                allGood = false;
            }
        }
        //coast should all be taken from the same parent
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
    else if (printMask){
        cout << endl;
    }

    crossOver_averageRatio(mask);
    //all the values in the indices in crossOver_average should contain threes
    for (int i = 0; i < OPTIM_VARS; i++){
        if (printMask){
            cout << mask[i] << " ";
        }
        if (mask[i] != 4){
            averageRatioGood = false;
            allGood = false;
        }
    }
    if (!averageRatioGood){
        cout << "crossOver_average generates the mask incorrectly" << endl;
    }
    else if (printMask){
        cout << endl;
    }

    return allGood; //returns whether or not there were errors deteced
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
    utcConstants->default_mutation_chance = 0.0; 

    for (int i = 0; i < 2; i++){
        childrenCreated = 0;

        //Create a mask to determine which of the child's parameters will be inherited from which parents
        std::vector<int> mask;
        for (int j = 0; j < OPTIM_VARS; j++) //initializes mask to average in case it gets stuck - also to get the correct length preset
        {
            mask.push_back(AVG); 
        }

        double annealing = 0.0; //we don't want any mutations yet

        //Generate a pair of children based on the random cross over mask method from a pair of parents
        //Generate the base mask
        crossOver_wholeRandom(mask, rng);

        //Generate a pair of children based on the mask
        generateChildrenPair(parent1, parent2, children, utcConstants->num_individuals, mask, annealing, rng, childrenCreated, gen, utcConstants);

        //the "magic" 1 represents the first test, crossover_wholeRandom
        if (!checkReasonability(children[childrenCreated -2], children[childrenCreated-1], mask, parentVals, 1, utcConstants, printThings)){
            noErrors = false;
        }

        //Generate a pair of children from a pair of parents based on the average mask method
        //Generate the base mask
        crossOver_average(mask);

        //Generate a pair of children based on the mask
        generateChildrenPair(parent1, parent2, children, utcConstants->num_individuals, mask, annealing, rng, childrenCreated, gen, utcConstants);

        //the "magic" 2 represents the second test, crossover_average
       if (!checkReasonability(children[childrenCreated -2], children[childrenCreated-1], mask, parentVals, 2, utcConstants, printThings)){
            noErrors = false;
        }
        //Generate a pair of children from a pair of parents based on the bundled random mask method
        //Generate the base mask
        crossOver_bundleVars(mask, rng);

        //Generate a pair of children based on the mask
        generateChildrenPair(parent1, parent2, children, utcConstants->num_individuals, mask, annealing, rng, childrenCreated, gen, utcConstants);
        
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
bool checkReasonability(const Child& c1, const Child& c2, std::vector<int> & mask, double* parentsValues, int whichMethod, cudaConstants * utcConstants, bool printThings){
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
                cout << "\tError with ";
                skipPrint = false; //whether or not valid messages are being printed, it needs to print error message
                noErrors = false;
            }
            else if (printThings){
                cout << "\tExpected values for ";   
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
                cout << "\tError with ";
                skipPrint = false; //whether or not valid messages are being printed, it needs to print error message
                noErrors = false;
            }
            else if (printThings){
                cout << "\tExpected values for ";   
            }
            else{
                skipPrint = true;
            }
        }
        //the only time it would not be an error for a child to have AVG as a value in its mask is when the method being used is crossOver_average (whichMethod == 2)
        else if (mask[i] == AVG_RATIO && whichMethod == 2){
            //if the value is supposed to be an average, ensures it actually is -> only the first child generated is a true average, the second is a weigthed average
            if (getParamStuff(i,c1) != (parentsValues[parValIndex]+parentsValues[parValIndex+1])/2){
                cout << "\tError with ";
                skipPrint = false; //whether or not valid messages are being printed, it needs to print error message
                noErrors = false;
            }
            //ensures that the random average ends up with something that is at least within the bounds of where it should be 
            else if ((getParamStuff(i,c2) > parentsValues[parValIndex] && getParamStuff(i,c2) > parentsValues[parValIndex+1]) || (getParamStuff(i,c2) < parentsValues[parValIndex] && getParamStuff(i,c2) < parentsValues[parValIndex+1])){
                cout << "\tError with ";
                skipPrint = false; //whether or not valid messages are being printed, it needs to print error message
                noErrors = false;
            }
            
            else if (printThings){
                cout << "\tExpected values for ";   
            }
            else{
                skipPrint = true;
            }
        }
        //if the values in the mask does not correspond to PARTNER1,PARTNER2, or AVG_RATIO then there is an issue with the mask
        else{
            cout << "\tError with mask for ";
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
double getParamStuff(const int & correspondingOffset, const Child& aChild){
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

void twentyAdultsPosAndSpeedDiffMade(bool printThings, std::vector<Adult>& allAdults, cudaConstants* utcConstants){
    utcConstants->num_individuals = 20;

    Child* genZero = new Child[utcConstants->num_individuals]; 
    allAdults.clear();

    //made 20 children with semi-random speed and position differences and tripTimes
    //the tripTime has no actual correlation with the posDiff and speedDiff
    //some of these tripTimes values are not correctly
    genZero[0] = Child(40000000.0, 0.02, 0.0043); //about 1.27 years  
    genZero[1] = Child(35000000.0, 0.048, 0.00234); //about 1.12 years 
    genZero[2] = Child(41000000.0, 0.299, 0.0034); //about 1.30 years
    genZero[3] = Child(32000000.0, 0.025, 0.00354); //about 1.01 years 
    genZero[4] = Child(47000000.0, 0.122, 0.0034); //about 1.49 years 
    genZero[5] = Child(43000000.0, 0.02, 0.0034); //about 1.36 years 
    genZero[6] = Child(37400000.0, 0.012, 0.0042); //about 1.18 years 
    genZero[7] = Child(38000000.0, 0.14, 0.00192); //about 1.20 years 
    genZero[8] = Child(44300000.0, 0.053, 0.00414); //about 1.40 years 
    genZero[9] = Child(45200000.0, 0.098, 0.00432); //about 1.43 years 
    genZero[10] = Child(40000010.0, 0.02, 0.0043); //about 1.27 years  
    genZero[11] = Child(35000020.0, 0.048, 0.00234); //about 1.12 years 
    genZero[12] = Child(41000030.0, 0.299, 0.0034); //about 1.30 years
    genZero[13] = Child(32000040.0, 0.025, 0.00354); //about 1.01 years 
    genZero[14] = Child(47000050.0, 0.122, 0.0034); //about 1.49 years 
    genZero[15] = Child(43000060.0, 0.02, 0.0034); //about 1.36 years 
    genZero[16] = Child(37400070.0, 0.017, 0.0042); //about 1.18 years 
    genZero[17] = Child(38000080.0, 0.14, 0.00192); //about 1.20 years 
    genZero[18] = Child(44300090.0, 0.02, 0.0034); //about 1.40 years 
    genZero[19] = Child(45200100.0, 0.299, 0.0034); //about 1.43 years

    //all the children are converted intp Adults
    for (int i = 0; i < utcConstants->num_individuals; i++){
        allAdults.push_back(Adult(genZero[i]));
    }

    //ranks and gives distances to all the adults and then sorts them
    giveRank(allAdults, utcConstants);
    giveDistance(allAdults, utcConstants);
    std::sort(allAdults.begin(), allAdults.end(), rankDistanceSort);

    //if printThings is true, it will print all the individuals in rankDistanceSort order
    if (printThings){
        for (int i = 0; i < utcConstants->num_individuals; i++){
            cout << "allAdults[" << i << "]: " << allAdults[i].unitTestingRankDistanceStatusPrint() << "- posDiff: " << allAdults[i].posDiff << ", speedDiff: " << allAdults[i].speedDiff << " & tripTime: " << allAdults[i].startParams.tripTime << endl;
        }
        cout << endl;
    }

    delete[] genZero; //cleans up dynamic memory
}

//checks to make sure children are being generated the way we think they are
bool verifyChildrenFromCrossover(std::mt19937_64& rng, bool printThings, cudaConstants* utcConstants){
    //holds the potential parents
    std::vector<Adult> oldAdults;
    //holds these potential parents based on whether or not they are duplicates or unique/original values
    std::vector<Adult> parents;

    //this vector is necessary because now findDuplicates takes in two vectors, so since we only want to find the duplicates...
    //  ...in oldAdults, we pass this empty vector in along with oldAdults
    std::vector<Adult> empty;

    Child occupiesLastSpace;

    //holds the status of the code - whether or not we ran into any errors -> this will be the output
    bool noErrors = true;

    //the power of two which is dividing num_individuals to get the survivpr count
    int exponent = 0;
    int numDupes = 0;
    int generation = 1;

    bool printParents = false; //set the print status to false because we do not need to print the 20 adults again before they are modified

    //sets the population size to twenty and fills old Adults with 20 Adults
    twentyAdultsPosAndSpeedDiffMade(printParents, oldAdults, utcConstants);
    
    //starts the number of survivors at the entire population of individuals
    utcConstants->survivor_count = utcConstants->num_individuals;

    //holds the children to be generated
    //this must be initialized AFTER num_individuals is updated to 20
    Child* children = new Child[utcConstants->num_individuals]; 

    //loops 4 times the number of survivors is n/(2^exponent) -> survivor_count = n (20), n/2 (10), n/4 (5), n/8 (2)
    while (exponent < 4){
        //if it going to be printing information, it tell the user the survivor size and population size 
        if (printThings){
            cout << "\nSurvivor population: " << utcConstants->survivor_count << " out of " << utcConstants->num_individuals << " potential parents" << endl;
            cout << "With a pool of entirely unique individuals..." << endl;
        }
        bool goodThisTime = true; //initially assumes the test will pass
        printParents = false; //set the print status to false because we do not need to print the 20 adults again before they are modified

        //sets the population size to twenty and fills old Adults with 20 Adults
        twentyAdultsPosAndSpeedDiffMade(printParents, oldAdults, utcConstants);

        //fills oldAdults with the 20 adults who only have posDiff, speedDiff, and tripTime uniquely assigned to them
        //twentyAdultsPosAndSpeedDiffMade(printParents, oldAdults, utcConstants);
        empty.clear();

        //assigns a duplicate status to all duplicates in oldAdults
        findDuplicates(oldAdults, empty, utcConstants, utcConstants->anneal_initial); 

        //a vector of unique posDiffs to replace the posDiff of any duplicate with
        std::vector<double> posDiffs;
        posDiffs.push_back(0.23);
        posDiffs.push_back(0.153);
        posDiffs.push_back(0.284);
        posDiffs.push_back(0.032);
        posDiffs.push_back(0.835);
        posDiffs.push_back(0.556);
        posDiffs.push_back(0.679);
        posDiffs.push_back(0.332);
        posDiffs.push_back(0.026);

        //a vector of unique speedDiffs to replace the speedDiff of any duplicate with
        std::vector<double> speedDiffs;
        speedDiffs.push_back(0.009221);
        speedDiffs.push_back(0.003865);
        speedDiffs.push_back(0.008965);
        speedDiffs.push_back(0.007342);
        speedDiffs.push_back(0.00856);
        speedDiffs.push_back(0.001356);
        speedDiffs.push_back(0.0014368);
        speedDiffs.push_back(0.008365);
        speedDiffs.push_back(0.009243);
        
        //make all the Adults in oldAdults unique individuals
        for (int i = 0; i < oldAdults.size(); i++){
            //if an individual is a duplicate, replace both its posDiff and speedDiff with different values to make it unique
            if (oldAdults[i].errorStatus == DUPLICATE){
                oldAdults[i].posDiff = posDiffs[posDiffs.size()-1];
                posDiffs.pop_back();
                oldAdults[i].speedDiff = speedDiffs[speedDiffs.size()-1];
                speedDiffs.pop_back();
            }
        }

        //ensure empty is actually empty, then tries to locate duplicates in oldAdukts
        //there should no longer be any duplicates left
        empty.clear();
        findDuplicates(oldAdults, empty, utcConstants, utcConstants->anneal_initial);

        //ensures all the adults in oldAdults are unique
        for (int i = 0; i < oldAdults.size(); i++){
            //if there is a duplicate, it prints an error message and will print oldAdults so the user can examine the issue
            if (oldAdults[i].errorStatus == DUPLICATE){
                cout << i << endl;
                if (speedDiffs.size() == 0 || posDiffs.size() == 0){
                    cout << "There were not enough unique speedDiffs and posDiffs to fix all the duplicate values or findDuplicates isn't working as expected" << endl;
                }
                else {
                    cout << "Other unknown error having to do with proper posDiffs and speedDiffs not being assigned" << endl;
                }
                printParents = true;
            }
        }

        //ensures parents is prepareed to take the new parents
        parents.clear();
        //separates oldAdults into parents (full of best N/4 individuals) 
        fillParents(oldAdults, parents, generation, utcConstants);

        for (int i = 0; i < utcConstants->num_individuals; i++){
            if (i < utcConstants->survivor_count){
                if(parents[i].errorStatus == DUPLICATE){
                    numDupes++;
                }
            }
            //prints all the identifying information (including the error status for each Adult in parents)
            if (printParents){
                cout << "oldAdults[" << i << "]: " << oldAdults[i].unitTestingRankDistanceStatusPrint() << "- posDiff: " << oldAdults[i].posDiff << ", speedDiff: " << oldAdults[i].speedDiff << " & tripTime: " << oldAdults[i].startParams.tripTime << endl;
            }
        }
        if (printParents){
            cout << endl;
            printParents = false; //it will not do this large print next time 
        }
        //no mutations should occur because anneal_initial is set to 0, as is the default_mutation_chance
        utcConstants->anneal_initial = 0;
        utcConstants->default_mutation_chance = 0;

        //set the generation at 1 because this is the first generation and that should have little influence on the generation
        //the number of children to be generated using this method is num_individuals because there should be no duplicates at this point
        generateChildrenFromCrossover(parents, children, utcConstants->num_individuals, rng, utcConstants->anneal_initial, generation, utcConstants);
        goodThisTime = cfcAnswersMatchExpectations(utcConstants->num_individuals, children, utcConstants, parents, printThings, numDupes);
        //if the function that is meant to check if the answers match our expectations returns false, there are errors
        if (!goodThisTime){
            //identifies that there were errors if printThings is true
            if (printThings){
                cout << "\n...generateChildrenFromCrossover() does not work as expected" << endl;
            }
            noErrors = false;
        }

        //resets after each time it is calculated 
        numDupes = 0;

        //identifies the next test being done
        if (printThings){
            cout << "In a pool that is a mix of unique individuals and duplicates..." << endl;
        }

        oldAdults.clear();
        //refills oldAdults with 11 unique individuals and 9 duplicates 
        twentyAdultsPosAndSpeedDiffMade(printParents, oldAdults, utcConstants);
        empty.clear();
        findDuplicates(oldAdults, empty, utcConstants, utcConstants->anneal_initial); 
        //ensures parents is prepareed to take the new parents
        parents.clear();
        //separates oldAdults into parents (full of best N/4 individuals) 
        fillParents(oldAdults, parents, generation, utcConstants);

        for (int i = 0; i < utcConstants->survivor_count; i++){
            if(parents[i].errorStatus == DUPLICATE){
                numDupes++;
            }
            //prints all the identifying information (including the error status for each Adult in parents)
            if (printParents){
                cout << "parents[" << i << "]: " << parents[i].unitTestingRankDistanceStatusPrint() << "- posDiff: " << parents[i].posDiff << ", speedDiff: " << parents[i].speedDiff << " & tripTime: " << parents[i].startParams.tripTime << endl;
            }
        }
        if (printParents){
            cout << endl;
            printParents = false; //it will not do this large print next time 
        }
        
        //generates individuals based on solely the original/unique individuals from the population
        //generation is set to two to distinguish individuals created in this second call from those created in the first
        generateChildrenFromCrossover(parents, children, utcConstants->num_individuals, rng, utcConstants->anneal_initial, generation, utcConstants);
        goodThisTime = cfcAnswersMatchExpectations(utcConstants->num_individuals, children, utcConstants, parents, printThings, numDupes);
        //if the function that is meant to check if the answers match our expectations returns false, there are errors
        if (!goodThisTime){
            //identifies that there were errors if printThings is true
            if (printThings){
                cout << "\n...generateChildrenFromCrossover() does not work as expected" << endl;
            }
            noErrors = false;
        }

        //resets this so it can accurately hold the number of duplicates
        numDupes = 0;

        // ===== NOTE =====
        // Attempted having a pool of entirely duplicates 
        // there was only 1 individual in parent and generateChildrenFromCrossover could not handle it
        // it began an infinite loop of printing "In generateChildrenFromCrossover() we've cycled through the parentPool twice"
        // the code for this test is commented out below
        // TODO: What should I do with this code? We should never get this situation, but I thought I should test it

        if (printThings){
            cout << "In a pool full of duplicates of one individual..." << endl;
        }

        //fill oldAdults entirely with duplicates 
        for (int i = 0; i < oldAdults.size(); i++){
            oldAdults[i] = Adult(Child(3.8e+07, 0.014, 0.00192));
        }

        //locates and separates the duplicates
        empty.clear();
        findDuplicates(oldAdults, empty, utcConstants, utcConstants->anneal_initial); 
        //ensures parents is prepareed to take the new parents
        parents.clear();
        //separates oldAdults into parents (full of best N/4 individuals) 
        fillParents(oldAdults, parents, generation, utcConstants);

        for (int i = 0; i < utcConstants->survivor_count; i++){
            if(parents[i].errorStatus == DUPLICATE){
                numDupes++;
            }
            //prints all the identifying information (including the error status for each Adult in parents)
            if (printParents){
                cout << "parents[" << i << "]: " << parents[i].unitTestingRankDistanceStatusPrint() << "- posDiff: " << parents[i].posDiff << ", speedDiff: " << parents[i].speedDiff << " & tripTime: " << parents[i].startParams.tripTime << endl;
            }
        }
        if (printParents){
            cout << endl;
            printParents = false; //it will not do this large print next time 
        }

        //generates individuals based on solely the original/unique individual from oldAdults
        //can the code handle only one parent? If so, what does it do?
        //generation is set to three to distinguish individuals created in this third call from those created in the first and second
        generateChildrenFromCrossover(parents, children, utcConstants->num_individuals, rng, utcConstants->anneal_initial, generation, utcConstants);
        goodThisTime = cfcAnswersMatchExpectations(utcConstants->num_individuals, children, utcConstants, parents, printThings, numDupes);

        if (!goodThisTime){
            noErrors = false;
        }

        numDupes = 0;
        utcConstants->survivor_count /= 2; //Halves the survivor population
        exponent++; //increases the exponent
    }
    delete[] children; //deallocates dynamic memory
    
    return noErrors; //returns whether or not the tests were successful
}

// Compares the results from verifyChildrenFromCrossover with the correct values
bool cfcAnswersMatchExpectations(const int & numChildren, const Child* childrenGenerated, const cudaConstants* utcConstants, std::vector<Adult> & parents, bool printThings, int duplicateCount){
    bool noErrors = true; //initially assumes that all the values passed in will meet our expectations
    bool skipPrint = !printThings;

    //if there numChildren does not fill the entire childrenGenerated array, ensures that it did not overwrite a value that was not supposed to change
    if (printThings){
        if (duplicateCount == (parents.size() - 1)){
            cout << "All the parents are duplicates of one individual. So, there will not be a print out as the results would not be very interesting " << endl;
            skipPrint = true;
        }
        if (parents.size() < 3){
            cout << "There are only two parents, so there will be little variation in the Children, so nothing will be printed " << endl;
            skipPrint = true;
        }
    }
    
    int currSet = 0; //the set of children this is examining
    while(currSet*crossoverChildrenCount < numChildren){
        bool correctParents = false; //true if we can be fairly sure a Child had two unique and valid parents
        bool ableToCheckBothParents = false; //this is only false if we have one more child than would make a complete set (7, 13, 19, etc)
        int prnt1 = -1, prnt2 = -1; //represents the index of a Child's parents
        //if a pair of parents have generated 3 or more children
        if (currSet*crossoverChildrenCount + fullAvgOffset < numChildren){ 
            ableToCheckBothParents = true; //we will be able to identify two unique parents
            for (int i = 0; i < utcConstants->survivor_count; i++){
                //a parent must have a partner to generate a Child, so it is counted as an error if a parent generates a Child by itself
                for (int j = i+1; j < utcConstants->survivor_count; j++){ 
                    if ((parents[i].startParams.tripTime + parents[j].startParams.tripTime)/2 == childrenGenerated[currSet*crossoverChildrenCount + fullAvgOffset].startParams.tripTime){
                        correctParents = true;
                        prnt1 = i;
                        prnt2 = j;
                        //TODO: Run the code and make sure the output of this is as expected/looks nice
                        //prints the tripTimes of the parents and of any children they've made
                        if (!skipPrint){
                            cout << std::fixed << std::setprecision(1) << "Parent 1's tripTime: " << parents[prnt1].startParams.tripTime << "\t\t Parent 2's tripTime: " << parents[prnt2].startParams.tripTime << endl;

                            int kids = 0; //how many children this parent set had

                            //loops until it has printed the tripTime for every Child a couple made
                            while(currSet*crossoverChildrenCount + kids < numChildren && kids < crossoverChildrenCount){
                                cout << std::fixed << std::setprecision(1) << "Child " << kids+1 << "'s tripTime: " << childrenGenerated[currSet*crossoverChildrenCount+kids].startParams.tripTime;
                                //prints two Children's tripTimes per line, so only ends the line every other time
                                if (kids%3 == 2){
                                    cout << endl;
                                }   
                                else {
                                    cout << "\t\t ";
                                    if (currSet*crossoverChildrenCount + kids+1 == numChildren){
                                        cout << endl;
                                    }
                                }
                                kids++; //increments kids and checks if the parents created another one
                            }  
                            cout << endl; //prints an extra blank line between every output
                        }

                        //sets i and j to their limits to force it to break from the loop
                        //otherwise when there are duplicates, it will print practically the same parent and child set multiple times
                        i = parents.size();
                        j = parents.size();
                    }
                }
            }
        }
        //if a pair of parents generated only 1 or 2 children
        else{
            //attempts to find a parent whose tripTime values match the first child's tripTimes values because this child is generated by the crossOver_wholeRandom method
            for (int i = 0; i < parents.size(); i++){
                //if a parent's tripTime matches the first Child's tripTime, it was a parent of this incomplete set of children
                if ((parents[i].startParams.tripTime == childrenGenerated[currSet*crossoverChildrenCount].startParams.tripTime) && prnt1 == -1){
                   prnt1 = i;
                }
                //if there is a second child that can be checked, we may be able to identify the other parent
                if (currSet*crossoverChildrenCount + 1 < numChildren){
                    ableToCheckBothParents = true; //there is a chance we can identify both parents
                    //if a parent's tripTime matches this second child's tripTime it is also assumed to be a parent of the incomplete set
                    if ((parents[i].startParams.tripTime == childrenGenerated[currSet*crossoverChildrenCount + 1].startParams.tripTime) && prnt1 != i){
                        prnt2 = i;
                    }
                }
            }
            //if printThings is true, it will print the parents' and their children's tripTimes
            if (!skipPrint){
                cout << std::fixed << std::setprecision(1) << "Parent 1's tripTime: " << parents[prnt1].startParams.tripTime << "\t\t Parent 2's tripTime: ";
                //if the code identified two distinct parnets, print each of their tripTimes and their Children's
                //to be able to identify two unique parents this way, there must be two Children so prints both
                if (ableToCheckBothParents && prnt1 != prnt2){
                    cout << std::fixed << std::setprecision(1) << parents[prnt2].startParams.tripTime << endl;
                    cout << std::fixed << std::setprecision(1) << "Child 1's tripTime: " << childrenGenerated[currSet*crossoverChildrenCount].startParams.tripTime;
                    cout << std::fixed << std::setprecision(1) << "\t\t Child 2's tripTime: " << childrenGenerated[currSet*crossoverChildrenCount +1].startParams.tripTime << endl;
                }
                //if it somehow fails to detect to distinct parents, but has two different Children (this should never happen)
                //prints both Children's tripTimes and identifies that the second parent's tripTime is unknown
                else if (ableToCheckBothParents){
                    cout << "UNKNOWN" << endl;
                    cout << std::fixed << std::setprecision(1) << "Child 1's tripTime: " << childrenGenerated[currSet*crossoverChildrenCount].startParams.tripTime;
                    cout << std::fixed << std::setprecision(1) << "\t\t Child 2's tripTime: " << childrenGenerated[currSet*crossoverChildrenCount +1].startParams.tripTime << endl;
                }
                //if there was only one child, it could not identify the second parent
                //prints a message that it could not identify the second parent's tripTime 
                //and prints the known Child's tripTime
                else {
                    cout << "UNKNOWN" << endl;
                    cout << std::fixed << std::setprecision(1) << "Child 1's tripTime: " << childrenGenerated[currSet*crossoverChildrenCount].startParams.tripTime << endl;
                }
                cout << endl; //prints an extra blank line at the end
            }
        }

        //if prnt1 was given a value
        if (prnt1 > -1){
            //it will identify whether it only checked one or both parents
            if (ableToCheckBothParents){
                //if parent 2 has been identified and is unique, it has found both parents
                if (prnt2 > -1 && prnt2 != prnt1){
                    correctParents = true;
                }
                //otherwise, there was an error
                else{
                    correctParents = false;
                }
            }
            //if only one parent was identified, it cannot conclusively say whether or not the set would have the correct parents...
            //  ...since tripTimes was the only variable that was set
            else{
                cout << "Unable to determine if children were correctly generated from parents" << endl;
                correctParents = true; 
            }
        }
        //if prnt1 is not even identified, it has not correctly identified the parents
        else {
            correctParents = false;
        }
        //if it did not correctly identify the parents, then there are errors
        if (!correctParents){
            if (printThings){
                cout << "Could not correctly identify the parent" << endl;
            }
            noErrors = false;
        }
        currSet++; //goes on to checking the next set of children
    }
    return noErrors; //if all children had their parents identified and they seem good, this returns true
}

bool testingCopyOfNewGen(bool printThings, std::mt19937_64 & rng, cudaConstants* utcConstants){
    //originally assumes that all the functions are working the way we expect them to
    bool allWorking = true;
    //starts the number of survivors and individuals at 20 so we can use all the duplicates to make children
    utcConstants->num_individuals = 20;
    utcConstants->survivor_count = 20;
    //ensures anneal_initial is set as well as the default and duplicate mutation chances
    utcConstants->anneal_initial = 0.01;
    utcConstants->default_mutation_chance = 0.0;

    //holds the potential parents
    std::vector<Adult> oldAdults;
    std::vector<Adult> newAdults;

    //holds the newAdults once they have been converted back into children for verification
    Child* newChildren = new Child[utcConstants->num_individuals]; 

    //sets the generation to 1 because generation should not have a major impact on any of these function
    int generation = 1;
    int duplicateCount = 0;

    cout << endl;

    //the hard-coded false is just supposed to stop this from printing all 20 Adults again
    twentyAdultsPosAndSpeedDiffMade(false, oldAdults, utcConstants);

    //UTCopyOfNewGen returns the number of individuals that should have been generated by crossOver
    int endSpotHere = UTCopyOfNewGeneration(oldAdults, newAdults, utcConstants->anneal_initial, generation, rng, utcConstants);

    //the two verification functions (cfcAnswersMatchExpectations and properMutations need an array of Children as their inputs)
    convertBackToChildren(newAdults, newChildren, utcConstants);

    //Holds the parents and duplicates -> these should be the same as the ones with the same name in UTCopyOfNewGeneration
    std::vector<Adult> parents; 
    parents.clear();
    //separates oldAdults into parents (full of best N/4 individuals) 
    fillParents(oldAdults, parents, generation, utcConstants);

    for (int i = 0; i < parents.size(); i++){
        if (parents[i].errorStatus == DUPLICATE){
            duplicateCount++;
        }
    }

    //if the children from crossover were not correctly generted, it makes the values that will be returned false 
    if (!cfcAnswersMatchExpectations(endSpotHere, newChildren, utcConstants, parents, printThings, duplicateCount)){
        allWorking = false;
    }

    //if both tests worked as expected, returns TRUE
    return allWorking;
}

//converts Adults to Children so they can be put into the existing verification algorithms
void convertBackToChildren(std::vector<Adult>& newAdult, Child* newChildren, const cudaConstants* utcConstants){
    //both newAdults and newChildren should be num_individuals size
    for (int i = 0; i < utcConstants->num_individuals; i++){
        //uses a unit testing Child constructor to create Children with pretty much all the same elements as the Adult has
        newChildren[i] = Child(newAdult[i].startParams, newAdult[i].posDiff, newAdult[i].speedDiff, newAdult[i].progress, newAdult[i].avgParentProgress, newAdult[i].birthday);
        //if the posDiffs are not the same that is a good sign that something really isn't working with the Child constructor
        if (newChildren[i].posDiff != newAdult[i].posDiff || newChildren[i].speedDiff != newAdult[i].speedDiff){
            cout << "This needs a better constructor " << endl;
        }
    }
}

//copied from genetic_algorithm.cpp on 6/29/22 at 5:59PM and then removed all the callRK and callRK dependent things
//made it return an int representing where the children from the unique individuals ended and the duplicates started
int UTCopyOfNewGeneration(std::vector<Adult> & oldAdults, std::vector<Adult> & newAdults, const double & annealing, const int & generation, std::mt19937_64 & rng, const cudaConstants* utcConstants){
    //Vector that will hold the adults who are potential parents
    //The criteria for being a parent is being in the top survivor_count number of adults in the oldAdult pool
    std::vector<Adult> parents; 

    parents.clear();

    //separates oldAdults into parents (full of best N/4 individuals) 
    fillParents(oldAdults, parents, generation, utcConstants);

    //Import number of new individuals the GPU needs to fill its threads
    //Create a newChildren function to fill in with children generated from oldAdults
    Child* newChildren = new Child[utcConstants->num_individuals]; 

    //uses parents to make children for the next generation
    makeChildren(parents, newChildren, annealing, generation, rng, utcConstants);
    
    //Now that the children have been simulated, convert the children into adults
    //First, it will calculate the right pos and speed diffs for the children
    //This will also put the converted children into the newAdults vector
    convertToAdults(newAdults, newChildren, utcConstants); 

    //Free the pointer's memory
    delete[] newChildren; 
    return (utcConstants->num_individuals);
}

//unit testing creating an entire generation of individuals
bool firstFullGen(std::mt19937_64& rng, cudaConstants * utcConstants, bool printThings){
    utcConstants->num_individuals = 10;
    utcConstants->survivor_count = 3;

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
    paramsForIndividuals[6] = rkParameters<double>(37400000.0, elems, coeffs); //about 1.18 years 
    paramsForIndividuals[7] = rkParameters<double>(38000000.0, elems, coeffs); //about 1.20 years 
    paramsForIndividuals[8] = rkParameters<double>(44300000.0, elems, coeffs); //about 1.40 years 
    paramsForIndividuals[9] = rkParameters<double>(45200000.0, elems, coeffs); //about 1.43 years 
    
    //turns the rkParameters into children to make up genZero
    for (int i = 0; i < utcConstants->num_individuals; i++){
        //all these children were created in the first generation and have no parents, so the parent 
        genZero[i] = Child(paramsForIndividuals[i], utcConstants, 1, 0);
    }

    std::vector<Adult> parents;
    firstGeneration(genZero, parents, utcConstants); //genZero and is turned into parents using firstGeneration which has been unit tested previously

    //the parents are sorted so they can be selected to be parents for a new generation
    giveRank(parents, utcConstants); 
    std::sort(parents.begin(), parents.end(), rankSort); 
    giveDistance(parents, utcConstants); 
    std::sort(parents.begin(), parents.end(), rankDistanceSort);

    std::vector<Adult> youngGen;

    //creates young a young generation with an annealing rate of 0.0 and the generation is set to 1 (the random 0.0 and 1)
    newGeneration(parents, youngGen, 0.0, 1, rng, utcConstants);

    //verifies the children generated are as expected
    noErrors = verifyFullGen(youngGen, parents, utcConstants, printThings);

    //if the two vectors do not have the same size, there is definitely an issue
    if (youngGen.size() != parents.size() || parents.size() != utcConstants->num_individuals){
        noErrors = false;
    }

    //deallocates dynamic memory
    delete[] genZero;
    delete[] paramsForIndividuals;

    //if the first test passed successfully, tries making a lot a children from a small number of survivors
    if (noErrors){
        cout << "Test successful - trying creating a large number number of children from a small number of adults " << endl;
        noErrors = makeManyChildren(rng, youngGen, parents, utcConstants, printThings);
    }

    //returns noErrors
    return noErrors;
}

// A function that is used to verify that firstFullGen is working correctly
// Made this earlier before we started treating duplicates differently and 
//   so long as duplicates are handled the same as regular numbers this test should work
bool verifyFullGen(std::vector<Adult>& youngGen, std::vector<Adult>& possParents, const cudaConstants * utcConstants, bool printThings){
    bool noErrors = true;

    std::vector<int> parentIndexSets; //vector that holds the indices of sets of parents

    //int setsOfChildren = (utcConstants->num_individuals)/childrenPerSet;
    int setNum = -1; //starts at negative 1 because incremented on the first run

    int indexVecSize = 0;
    double tripTimeTol = 10e-6;

    //it checks which parents are being used and prints things if printThings is set to true
    for (int i = 0; i < utcConstants->num_individuals; i++){
        //when it has gone through a whole set generated by one set of parents, it increase setNum
        if (i%crossoverChildrenCount == 0){ 
            setNum++;
        }
        if (printThings){ //prints all the parents and children at each index
            if (i < possParents.size()){
                if (i < utcConstants->survivor_count){ //puts *| |* around any parent that has been selected to have children
                    cout << "*|";
                }
                cout << "parents[" << i << "]'s tripTime: ";
                cout << possParents[i].startParams.tripTime;
                if (i < utcConstants->survivor_count){
                    cout << "|*";
                }
                // else {
                //     cout << "\t";
                // }
            }
            else if (i > utcConstants->survivor_count){
                cout << "                               ";
            }   
            cout << "\t\t";

            cout << "youngGen[" << i << "]'s tripTime: ";
            cout << youngGen[i].startParams.tripTime << endl;
        }
        //the averaging function is the best way to determine which two parents created a child
        //none of the survivor count individuals should have the same average
        if (i == fullAvgOffset + crossoverChildrenCount*setNum){
            indexVecSize += 2; //each time this occurs, two parents should be added to the parentIndexSets vector
            //if the average tripTime is equal to the copied value of parent, just one parent created the child and it is a clone
            if (youngGen[i].startParams.tripTime == youngGen[6*setNum].startParams.tripTime){ 
                cout << "The parent did not have a partner and probably created clones" << endl;
                noErrors = false; //this is an error if a parent clones itself
                for (int j = 0; j < utcConstants->survivor_count; j++){
                    if (youngGen[i].startParams.tripTime == possParents[j].startParams.tripTime){
                        //cout << "Set " << setNum << " created by parent[" << j << "] and parent[" << j << "] " << endl;
                        parentIndexSets.push_back(j);
                        parentIndexSets.push_back(j);
                    }
                }
            }
            else{
                //otherwise, it searches different pairs of potentially chosen parents to determine which parents produced this offspring
                for (int j = 0; j < utcConstants->survivor_count; j++){
                    for (int k = j+1; k < utcConstants->survivor_count; k++){
                        //checks the two parents
                        if (youngGen[i].startParams.tripTime <= ((possParents[j].startParams.tripTime + possParents[k].startParams.tripTime)/2.0 + tripTimeTol) && youngGen[i].startParams.tripTime >= ((possParents[j].startParams.tripTime + possParents[k].startParams.tripTime)/2.0 - tripTimeTol)){
                            //cout << "Set " << setNum << " created by parent[" << j << "] and parent[" << k << "] " << endl;
                            parentIndexSets.push_back(j);
                            parentIndexSets.push_back(k);
                        }
                    }
                }
            }
            //if the parent has not been added, this is an issue and it prints error messages, 
            //  then pushes back -1 so parentIndexSets will be the proper size so we can detect if this happens again
            if(indexVecSize != parentIndexSets.size()){
                noErrors = false;
                //cout << indexVecSize << " vs " << parentIndexSets.size() << endl;
                cout << "Unknown parents for set #" << setNum << " of children" << endl;
                parentIndexSets.push_back(-1);
                parentIndexSets.push_back(-1);

            }
        }
        //if there are no other ways to determine the parents, determines t
        else if ((i == crossoverChildrenCount*setNum || i == 1 + crossoverChildrenCount*setNum) && crossoverChildrenCount*setNum + fullAvgOffset > youngGen.size()) {
            indexVecSize++;
            //attempts to find a parent whose tripTime values match the first child's tripTimes values because this child is generated by the crossOver_wholeRandom method
            for (int j = 0; j < utcConstants->survivor_count; j++){
                //if a parent's tripTime matches the first Child's tripTime, it was a parent of this incomplete set of children
                if (youngGen[i].startParams.tripTime == possParents[j].startParams.tripTime){
                    parentIndexSets.push_back(j);
                }
            }
            //if the parent has not been added, this is an issue and it prints error messages, 
            //  then pushes back -1 so parentIndexSets will be the proper size so we can detect if this happens again
            if(indexVecSize != parentIndexSets.size()){
                noErrors = false;
                //cout << indexVecSize << " vs " << parentIndexSets.size() << endl;
                cout << "Unknown parents for set #" << setNum << " of children" << endl;
                parentIndexSets.push_back(-1);

            }
        }
    }
    return noErrors;
}

//tries making children from just a couple of survivors - not currently working as far as I know...
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
    else if (printThings){
        cout << "Generated 65 individuals from 7 parents" << endl;
    }

    utcConstants->survivor_count = 5;
    newGeneration(possParents, youngGen, 0.0, 1, rng, utcConstants);
    //verifies 
    noErrors = verifyFullGen(youngGen, possParents, utcConstants, false);

    if (!noErrors){
        cout << "Problem with 65 individuals generated from 5 parents" << endl;
    }
    else if (printThings){
        cout << "Generated 65 individuals from 5 parents" << endl;
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



// bool runGeneticsUnitTests(bool printThings){
// // SETTING UP CUDA CONSTANTS TO BE USED BY OTHER FUNCTIONS 
//     //making cudaConstants to control what is going on in genetic algorithm while still using the original functions
//     cudaConstants* utcConstants = new cudaConstants(); 

//     //sets these tolerances so that the domination in stolenGiveRank works properly
//     //if these are not set, it gives the wrong rank to one Adult from each of the lower ranks (it promotes them to the next rank)
//     utcConstants->posDominationTolerance = 1.0e-14;
//     utcConstants->speedDominationTolerance = 1.0e-14;

//     //chose pos_threshold and speed_threshold so anything past the first significant figure in the clone separation test is trivial
//     utcConstants->pos_threshold = 1.0e-11;
//     utcConstants->speed_threshold = 1.0e-11;

//     // Seed used for randomization rng things, using seed 0 for consistancy / tracability 
//     utcConstants->time_seed = 0; 

//     //values taken directly from the genetic config file
//     utcConstants->anneal_initial = 0.0; // initial value for annealing, was 0.1 in config, set to 0 because no mutations wanted at first
//     utcConstants->anneal_final = 1.0e-7;   // final value for annealing, was 1.0e-7 in config, set to 0 because no mutations wanted at first
//     utcConstants->anneal_factor = 0.75;  // factor by which annealing is multiplied with when there is no change in the best individual over 100 generations -should be tested as-is

//     // The percentage for probability of mutating a gene in a new individual, called iteratively to mutate more genes until the check fails
//     // Starting this off at 0 because to ensure that there are no mutations initially to verify that children are generated as expected
//     // Changes to 1.0 later in the code so that we can ensure mutation is working
//     utcConstants->default_mutation_chance = 0; 
//     utcConstants->mutation_amplitude = 1;
   
//     // Used in mutate(), affects the scale of change for the respective parameter values, in conjunction with annealing
//     // Represents max bounds of mutation, mutation will never be +/- this value
//     utcConstants->gamma_mutate_scale = 3.14159; 
//     utcConstants->tau_mutate_scale = 1.570795; 
//     utcConstants->coast_mutate_scale = 3.14159;
//     utcConstants->triptime_mutate_scale = 1.0* SECONDS_IN_YEAR;
//     utcConstants->zeta_mutate_scale = 1.570795;
//     utcConstants->beta_mutate_scale = 1.570795;
//     utcConstants->alpha_mutate_scale = 3.14159;

//     // 0 is for no thruster, 1 is for NEXT ion thruster
//     // starts with no thruster to test that, then will be changed to have a thruster later on when necessary
//     utcConstants->thruster_type = 0; 

//     // Number of individuals in the pool -> chose to make this 10 so a generation size is managable
//     // Additionally a size of 10 shows what happens if a the number of children generated is not a factor of the needed generation size
//     // (E.g. 8 does not evenly divide 10 -> two pairings of parents will each produce 8 kids, which is 16 total and 16 > 10)  
//     utcConstants->num_individuals = 10; 

//     // Number of survivors selected, every pair of survivors creates 8 new individuals 
//     // Chose 3 because that is approximately a quarter of 10
//     // This also allows us to see if the shuffling and new parent pairing works as expected
//     // This will likely be changed in other sections of the code as well
//     utcConstants->survivor_count = 3;  

//     //Flag for what type of landing you want. Current modes: "soft"=2, "hard"=1
//     //testing for soft because as of 2022 that is the main mission type we are exploring 
//     //additionally, mission type should have little impact on the genetics algorithms
//     utcConstants->missionType = 2; 

//     //the maximum age an individual in this unit test can hit 
//     //chose 3 because this will allow us to minimize run while ensuring old Adults are properly removed
//     utcConstants->max_age = 3;

//     //STUFF SO RUNGE_KUTTACUDA.CU CAN RUN
//     utcConstants->thread_block_size =32;
//     utcConstants->orbitalPeriod = 3.772645011085093e+07; // orbital period time of 1999RQ36 Bennu (s) 

//     //ALL THE STUFF SO EARTH'S INITIAL POSITION CAN BE CALCULATED
//     utcConstants->triptime_max=2.0* SECONDS_IN_YEAR;
//     utcConstants->triptime_min=1.0* SECONDS_IN_YEAR;
//     utcConstants->timeRes=3600; // Earth Calculations Time Resolution Value
//     utcConstants->v_escape = 2162.4/AU; //sqrt of DART Mission c3energy because there are no c3energy values recorded for the bennu mission
//     //Bennu config values
//     utcConstants->r_fin_earth=9.857045197029908E-01;
//     utcConstants->theta_fin_earth=1.242975503287042;
//     utcConstants->z_fin_earth=-4.332189909686674E-05;
//     utcConstants->vr_fin_earth=-1.755004992027024E-09;
//     utcConstants->vtheta_fin_earth=2.019592304815492E-07;
//     utcConstants->vz_fin_earth=-2.312712519131594E-12;

//     // Various values that impact runge kutta
//     utcConstants->rk_tol=1e-12;
//     utcConstants->doublePrecThresh=1e-12;
//     utcConstants->max_numsteps=1000;
//     utcConstants->min_numsteps=400;

//     //New ga_crossover stuff
//     utcConstants->distanceTolerance=1.0e-14;
//     //default_mutation_factor=1.0         //Set at one to allow for the mutation scales below to be applied normally/to children at face value
//     //duplicate_mutation_factor=1.0;       //Set at 10; when duplicate adults are mutated, they will undergo stronger mutations


// // SETTING A RANDOM NUMBER GENERATOR (rng) TO BE USED BY FUNCTIONS
//     // This rng object is used for generating all random numbers in the genetic algorithm, passed in to functions that need it
//     std::mt19937_64 rng(utcConstants->time_seed);

//     //VERY complicated part of the code with some possibility of errors -> just needed for the child constructor with rkParameters and cConstants as its arguments
//     launchCon = new EarthInfo(utcConstants); 

// // CALLING THE DIFFERENT UNIT TESTING ALGORITHMS
//     bool allWorking = true;

//     std::vector<Adult> aAdults;
    
//     //firstParentsTest takes in the cuda constants and will verify that children can be converted into parents and sorted using rankDistanceSort
//     if (firstParentsTest(utcConstants, printThings)){
//         cout << "PASSED: Children can be converted to adults that can be sorted" << endl;
//     }
//     else{
//         cout << "FAILED: Children cannot be converted to adults that can be sorted" << endl;
//         allWorking = false;
//     }
//     //createMasks uses rng to create masks and ensure that their elements are as expected -> there are not AVGs where there should not be, etc
//     if (createMasks(rng, printThings)){
//         cout << "PASSED: All three masks were successfully generated as expected" << endl;
//     }
//     else{
//         cout << "FAILED: Not all the masks were successfully generated" << endl;
//         allWorking = false;
//     }
//     //uses rng and cuda constants to create children using different masks
//     //these children are converted to Adults and then tested to ensure their values match the expected values
//     if (makeChildrenWithDifferentMethods(rng, utcConstants, printThings)){
//         cout << "PASSED: Successfully made children using the three different methods" << endl;
//     }
//     else{
//         cout << "FAILED: Could not successfully make children using the three different methods" << endl;
//         allWorking = false;
//     }
//     if (verifyProperCloneSeparation(printThings, utcConstants)){
//         cout << "PASSED: Separating clones from unique numbers is working as expected" << endl;
//     }
//     else{
//         cout << "FAILED: Separating clones from unique numbers is not working as expected" << endl;
//         allWorking = false;
//     }
//     if (verifyChildrenFromCrossover(rng, printThings, utcConstants)){
//         cout << "PASSED: Children were generated successfully using the childrenFromCrossover function" << endl;
//     }
//     else{
//         cout << "FAILED: Children were not correctly generated using the childrenFromCrossover function" << endl;
//         allWorking = false;
//     }
//     //tested the compare functions within elements and rkParameters to see if they would be a good way to determine if mutations were working
//     //ultimately it turned out these functions have trouble detecting small changes in decimals 
//     //so even after fixing the issues they were having with nans, these functions did not end up being used in verifying the solutions from verifyChildrenFromMutation
//     if (testingCompareFunctions(printThings)){ 
//         cout << "PASSED: Compare operations for elements and rkParameters should work for checking if mutations work or not" << endl;
//     }
//     else{
//         cout << "FAILED: Do not use the compare operations for elements and rkParameters to check if mutations work or not" << endl;
//     }
//     if (verifyChildrenFromMutation(rng, printThings, utcConstants)){
//         cout << "PASSED: Children were generated successfully using the childrenFromMutation function" << endl;
//     }
//     else{
//         cout << "FAILED: Children were not correctly generated using the childrenFromMutation function" << endl;
//         allWorking = false;
//     }
//     if (testingCopyOfNewGen(printThings, rng, utcConstants)){
//         cout << "PASSED: newGeneration (without callRK) seems to be working as expected" << endl;
//     }
//     else{
//         cout << "FAILED: newGeneration (without callRK) does not seem to be working" << endl;
//     }
    
//     //creates a generation of parents and then creates children from these parents
//     //then these children are sent through a function that verifies their tripTime values 
//     //tripTime was chosen because it was the paramerter set in creating the parents and it is easiest to see if the children have the correct values for tripTime 
//     if (firstFullGen(rng, utcConstants, printThings)){
//         cout << "PASSED: Successfully made the first generation of adults from another set of adults" << endl;
//     }
//     else{
//         cout << "FAILED: Could not successfully make the first generation of adults from another set of adults" << endl;
//         allWorking = false;
//     }

//     delete utcConstants;
//     delete launchCon;
//     return allWorking;
// }

// //returns true if the first generation is generated
// bool firstParentsTest(const cudaConstants * utcConstants, bool printThings){
//     std::vector<Adult> parents;
//     Child* theChildren = new Child[utcConstants->num_individuals];

//     //Made 10 children and set them with semi-random position difference and speed difference values
//     //        that can easily be determined are either true or false
//     //These numbers are semi-random multiples of 10 and were entered in no specific order at this point
//     //This sort of simulates when the initial generation of oldAdults are created using random parameters 
//     //(but these don't change from run to run and are easier to do compare)
//     theChildren[0] = Child(1.5000, .2000); //
//     theChildren[1] = Child (1.2000, .9000); //
//     theChildren[2] = Child(1.8000, .3000); //
//     theChildren[3] = Child (.2000, .7000); //
//     theChildren[4] = Child(1.1000, .4000); //
//     theChildren[5] = Child (.3000, .9000); //
//     theChildren[6] = Child(.5000, 1.2000); //
//     theChildren[7] = Child (2.2000, 9.7000);
//     theChildren[8] = Child(.2000, 1.2000); //
//     theChildren[9] = Child (3.4000, .9000); //
//     for (int i = 0; i < utcConstants->num_individuals; i++){
//         parents.push_back(Adult(theChildren[i]));
//     }
//     //Ranks based on my hand calculation of ranks based on the version of giveRank from 6/16/22 
//     //Dominations are based on the dominations test
//     //Below is a list of the parameters of the above children and which rank they should belong to 
//     //    based on their position and speed differences
//     //The work is shown in a PDF on Teams in the 2022 folder
//     //Rank 1: (150, 20) distance- 1.176471; (20, 70) distance- 1.051546; (110, 40) distance - 0.247119
//     //Rank 2:  (180, 30) distance - 0.226501; (30, 90) distance - 0.108854; (20, 120) distance - 0.060340 
//     //Rank 3: (120, 90) distance - 1.111583; (50,120) distance - 0.117647
//     //Rank 4: (220,970) distance - 1.470588; (340,90) distance - 1.030928

//     giveRank(parents, utcConstants);
//     giveDistance(parents, utcConstants);
//     std::sort(parents.begin(), parents.end(), rankDistanceSort);
    
//     //The final result should be... (see the PDF for the work done to calculate this)
//     //Rank 1: (150, 20); (20, 70); (110, 40)
//     //Rank 2:  (180, 30); (30, 90); (20, 120)  
//     //Rank 3: (50,120); (120, 90)
//     //Rank 4: (220,970); (340,90)
    
//     if (printThings){
//         for (int i = 0; i < parents.size(); i++){
//             cout << "(" << parents[i].posDiff << "," << parents[i].speedDiff <<"): " <<parents[i].unitTestingRankDistanceStatusPrint() << " | ";
//         }
//         cout << endl;
//     }

//     delete[] theChildren;

//     //Ensures the size is correct and that the parents match their expected values
//     if (utcConstants->num_individuals == parents.size() && checkParentsTest(parents)){
//         return true;
//     }
//     else{
//         cout << "Could not turn children into parents to create a parent generation" << endl;
//         return false;
//     }

// }

// //A function that hold the expected values for parents test and that will compare the results
// //      from firstParentsTest to the correct order of the values 
// bool checkParentsTest(std::vector<Adult>& theResults){
//     std::vector<Adult> theAnswers;
//     theAnswers.push_back(Adult(Child(1.50, .20)));
//     theAnswers.push_back(Adult(Child(.20, .70)));
//     theAnswers.push_back(Adult(Child(1.10, .40)));
//     theAnswers.push_back(Adult(Child(1.80, .30)));
//     theAnswers.push_back(Adult(Child(.30, .90)));
//     theAnswers.push_back(Adult(Child(.20, 1.20)));
//     theAnswers.push_back(Adult(Child(.50,1.20)));
//     theAnswers.push_back(Adult(Child(1.20, .90)));    
//     theAnswers.push_back(Adult(Child(2.20,9.70)));
//     theAnswers.push_back(Adult(Child(3.40,.90)));

//     //if the two vectors are not the same size, we have a problem
//     if (theAnswers.size() != theResults.size()){
//         return false;
//     }

//     //loops through the results and compares its values to the expected values 
//     for (int i = 0; i < theResults.size(); i++){
//         //if either the posDiff or speedDiff is wrong, then the vector has been sorted wrong and it prints what the issie is
//         if (theResults[i].posDiff != theAnswers[i].posDiff || theResults[i].speedDiff != theAnswers[i].speedDiff){
//             cout << "The Adult here is (" << theResults[i].posDiff << "," << theResults[i].speedDiff << "), but it should be (" << theAnswers[i].posDiff << "," << theAnswers[i].speedDiff << ")" << endl;
//             return false;
//         }
//     }
//     return true;
// }

// //creates masks to ensure that they are generated properly
// bool createMasks(std::mt19937_64& rng, bool printMask){
//     bool wholeRandGood = true, averageGood = true, bundleVarsGood = true, averageRatioGood = true, allGood = true;
//     //Create a mask to determine which of the child's parameters will be inherited from which parents
//     //This is done the same way it is in newGeneration in ga_crossover.cpp
//     std::vector<int> mask;
//     for (int j = 0; j < OPTIM_VARS; j++) //initializes mask to average in case it gets stuck - also to get the correct length preset
//     {
//         mask.push_back(AVG); 
//     }

//     crossOver_wholeRandom(mask, rng);
//     //everything in crossOver_wholeRandom should be either 1s or 2s
//     for (int i = 0; i < OPTIM_VARS; i++){
//         if (printMask){
//             cout << mask[i] << " ";
//         }
//         if (mask[i] != 1 && mask[i] != 2){
//             wholeRandGood = false;
//             allGood = false;
//         }
//     }
//     if (!wholeRandGood){
//         cout << "\ncrossOver_wholeRandom generates the mask incorrectly" << endl;
//     }
//     else if (printMask){
//         cout << endl;
//     }

//     crossOver_average(mask);
//     //all the values in the indices in crossOver_average should contain threes
//     for (int i = 0; i < OPTIM_VARS; i++){
//         if (printMask){
//             cout << mask[i] << " ";
//         }
//         if (mask[i] != 3){
//             averageGood = false;
//             allGood = false;
//         }
//     }
//     if (!averageGood){
//         cout << "crossOver_average generates the mask incorrectly" << endl;
//     }
//     else if (printMask){
//         cout << endl;
//     }

//     crossOver_bundleVars(mask, rng);
//     //all the values in the indices in crossOver_bundleVars should contain sets of 1s or 2s
//     for (int i = 0; i < OPTIM_VARS; i++){
//         if (printMask){
//             cout << mask[i] << " ";
//         }
//         //if it is neither a 1 or a 2, there is definitely an error
//         if (mask[i] != 1 && mask[i] != 2){
//             bundleVarsGood = false;
//             allGood = false;
//         }
//         //gamma should all be taken from the same parent
//         if(i > GAMMA_OFFSET && i < GAMMA_OFFSET + GAMMA_ARRAY_SIZE){
//             if (mask[i] != mask[i-1]){
//                 cout << "\n" << i << ": " << mask[i] << " vs " << i-1 << ": " << mask[i-1] << endl;
//                 bundleVarsGood = false;
//                 allGood = false;
//             }
//         }
//         //tau should all be taken from the same parent
//         else if(i > TAU_OFFSET && i < TAU_OFFSET + TAU_ARRAY_SIZE){
//             if (mask[i] != mask[i-1]){
//                 cout << "\n" << i << ": " << mask[i] << " vs " << i-1 << ": " << mask[i-1] << endl;
//                 bundleVarsGood = false;
//                 allGood = false;
//             }
//         }
//         //coast should all be taken from the same parent
//         else if(i > COAST_OFFSET && i < COAST_OFFSET + COAST_ARRAY_SIZE){
//             if (mask[i] != mask[i-1]){
//                 cout << "\n" << i << ": " << mask[i] << " vs " << i-1 << ": " << mask[i-1] << endl;
//                 bundleVarsGood = false;
//                 allGood = false;
//             }
//         }
//     }
//     if(!bundleVarsGood){
//         cout << "crossOver_bundleVars generates the mask incorrectly" << endl;
//     }
//     else if (printMask){
//         cout << endl;
//     }

//     crossOver_averageRatio(mask);
//     //all the values in the indices in crossOver_average should contain threes
//     for (int i = 0; i < OPTIM_VARS; i++){
//         if (printMask){
//             cout << mask[i] << " ";
//         }
//         if (mask[i] != 4){
//             averageRatioGood = false;
//             allGood = false;
//         }
//     }
//     if (!averageRatioGood){
//         cout << "crossOver_average generates the mask incorrectly" << endl;
//     }
//     else if (printMask){
//         cout << endl;
//     }

//     return allGood; //returns whether or not there were errors deteced
// }

// bool makeChildrenWithDifferentMethods(std::mt19937_64& rng, cudaConstants * utcConstants, bool printThings){
//     const int expectedNumChildren = 6; //pairs of children should be generated 3 times, so the expected number of children that should be created are 6
//     int childrenCreated = 0;
//     //selected values for alpha (indices 0 & 1), beta (2 & 3), zeta (4 & 5), and tripTime that are within the acceptable range for rkParameters
//     //these values are near the edge of the acceptable range of values to make them easier to calculate
//     double parentVals[8] = {3.1414, -3.1414, 3.1414, 0.0, 1.5707, -1.5707,45000000.0, 35000000.0};
//     rkParameters<double> p1(parentVals[6], parentVals[0], parentVals[2], parentVals[4]); //using unit testing constructor with only tripTime, alpha, beta, and zeta
//     rkParameters<double> p2(parentVals[7], parentVals[1], parentVals[3], parentVals[5]);
//     Child par1(p1);
//     Child par2(p2);
//     Adult parent1(par1);
//     Adult parent2(par2);

//     Child* children = new Child[expectedNumChildren];

//     //sets the generation to 1 because generateChildPair takes in a generation
//     // 1 was chosen because this is the first generation and the generation number should not have a major impact on the code
//     int gen = 1; 

//     //while there have not been any errors detected in the codel this is set to true
//     bool noErrors = true;

//     //these individuals should not be mutated or else we cannot verify that generating individuals is working as expected
//     utcConstants->default_mutation_chance = 0.0; 

//     for (int i = 0; i < 2; i++){
//         childrenCreated = 0;

//         //Create a mask to determine which of the child's parameters will be inherited from which parents
//         std::vector<int> mask;
//         for (int j = 0; j < OPTIM_VARS; j++) //initializes mask to average in case it gets stuck - also to get the correct length preset
//         {
//             mask.push_back(AVG); 
//         }

//         double annealing = 0.0; //we don't want any mutations yet

//         //Generate a pair of children based on the random cross over mask method from a pair of parents
//         //Generate the base mask
//         crossOver_wholeRandom(mask, rng);

//         //Generate a pair of children based on the mask
//         generateChildrenPair(parent1, parent2, children, utcConstants->num_individuals, mask, annealing, rng, childrenCreated, gen, utcConstants);

//         //the "magic" 1 represents the first test, crossover_wholeRandom
//         if (!checkReasonability(children[childrenCreated -2], children[childrenCreated-1], mask, parentVals, 1, utcConstants, printThings)){
//             noErrors = false;
//         }

//         //Generate a pair of children from a pair of parents based on the average mask method
//         //Generate the base mask
//         crossOver_average(mask);

//         //Generate a pair of children based on the mask
//         generateChildrenPair(parent1, parent2, children, utcConstants->num_individuals, mask, annealing, rng, childrenCreated, gen, utcConstants);

//         //the "magic" 2 represents the second test, crossover_average
//        if (!checkReasonability(children[childrenCreated -2], children[childrenCreated-1], mask, parentVals, 2, utcConstants, printThings)){
//             noErrors = false;
//         }
//         //Generate a pair of children from a pair of parents based on the bundled random mask method
//         //Generate the base mask
//         crossOver_bundleVars(mask, rng);

//         //Generate a pair of children based on the mask
//         generateChildrenPair(parent1, parent2, children, utcConstants->num_individuals, mask, annealing, rng, childrenCreated, gen, utcConstants);
        
//         //the "magic" 3 represents the third test, crossover_bundleVars
//         if (!checkReasonability(children[childrenCreated -2], children[childrenCreated-1], mask, parentVals, 3, utcConstants, printThings)){
//             noErrors = false;
//         }

//         //if there were errors or not enough children were created on the first run with no thruster, won't bother testing with thruster
//         //the thruster should not make any difference in the parameters being checked by checkReasonability
//         if((childrenCreated != expectedNumChildren && i == 0) || !noErrors){
//             cout << "Could not even generate children correctly with no thruster" << endl;
//             delete[] children;
//             return false;
//         }
        
//         //once it has gone through once with no thruster, it turns the thuster on to make sure it works with the thruster too
//         utcConstants->thruster_type = 1; 
//     }

//     //if it executes correctly, deallocates dynamic memory and returns true
//     if(childrenCreated == expectedNumChildren && noErrors){
//         delete[] children;
//         return true;
//     }
//     //if there were errors, it deallocates dynamic memory and returns false
//     else{
//         cout << "Could not generate children correctly when there was thrust" << endl;
//         delete[] children;
//         return false;
//     }

// }

// //1 stands for crossOver_wholeRandom, 2 stands for crossOver_average, 3 stands for crossOver_bundleVars
// bool checkReasonability(const Child& c1, const Child& c2, std::vector<int> & mask, double* parentsValues, int whichMethod, cudaConstants * utcConstants, bool printThings){
//     bool noErrors = true; //initially there are no error detected yet
//     bool notFirstSwitch = false; //this is just here so the output turns out nicely
//     bool skipPrint = !printThings; //skipping printing is the opposite of printing
//     int parValIndex = 0;
//     for (int i = ALPHA_OFFSET; i <= TRIPTIME_OFFSET; i++){
//         //the mask was flipped, so compare the second child to it first, and the other child should have the other parent's parameters
//         //because flipping a mask sets PARTNER1 values to PARTNER2 and vice versa
//         if (mask[i] == PARTNER1 && whichMethod != 2){ 
//             //if the mask at this index corresponds to PARTNER1, the second child generated should have parent 1's value for this parameter
//             //and the first child should have the values the other parent (PARTNER2) has for this parameter
//             //if either of these are not the case, there is an issue
//             if (getParamStuff(i,c2) != parentsValues[parValIndex] || getParamStuff(i,c1) != parentsValues[parValIndex+1]){ 
//                 cout << parentsValues[parValIndex] << " should equal c2 " << getParamStuff(i,c2) << endl;
//                 cout << "\tError with ";
//                 skipPrint = false; //whether or not valid messages are being printed, it needs to print error message
//                 noErrors = false;
//             }
//             else if (printThings){
//                 cout << "\tExpected values for ";   
//             }
//             else{
//                 skipPrint = true;
//             }
//         }
//         else if (mask[i] == PARTNER2 && whichMethod != 2){
//             //if the mask at this index corresponds to PARTNER2, the second child generated should have parent 2's value for this parameter
//             //and the first child should have the values the other parent (PARTNER1) has for this parameter
//             //if either of these are not the case, there is an issue
//             if (getParamStuff(i,c2) != parentsValues[parValIndex+1] || getParamStuff(i,c1) != parentsValues[parValIndex]){
//                 cout << parentsValues[parValIndex+1] << " should equal c2 " << getParamStuff(i,c2) << endl;
//                 cout << "\tError with ";
//                 skipPrint = false; //whether or not valid messages are being printed, it needs to print error message
//                 noErrors = false;
//             }
//             else if (printThings){
//                 cout << "\tExpected values for ";   
//             }
//             else{
//                 skipPrint = true;
//             }
//         }
//         //the only time it would not be an error for a child to have AVG as a value in its mask is when the method being used is crossOver_average (whichMethod == 2)
//         else if (mask[i] == AVG_RATIO && whichMethod == 2){
//             //if the value is supposed to be an average, ensures it actually is -> only the first child generated is a true average, the second is a weigthed average
//             if (getParamStuff(i,c1) != (parentsValues[parValIndex]+parentsValues[parValIndex+1])/2){
//                 cout << "\tError with ";
//                 skipPrint = false; //whether or not valid messages are being printed, it needs to print error message
//                 noErrors = false;
//             }
//             //ensures that the random average ends up with something that is at least within the bounds of where it should be 
//             else if ((getParamStuff(i,c2) > parentsValues[parValIndex] && getParamStuff(i,c2) > parentsValues[parValIndex+1]) || (getParamStuff(i,c2) < parentsValues[parValIndex] && getParamStuff(i,c2) < parentsValues[parValIndex+1])){
//                 cout << "\tError with ";
//                 skipPrint = false; //whether or not valid messages are being printed, it needs to print error message
//                 noErrors = false;
//             }
            
//             else if (printThings){
//                 cout << "\tExpected values for ";   
//             }
//             else{
//                 skipPrint = true;
//             }
//         }
//         //if the values in the mask does not correspond to PARTNER1,PARTNER2, or AVG_RATIO then there is an issue with the mask
//         else{
//             cout << "\tError with mask for ";
//             noErrors = false;
//             skipPrint = false; //whether or not valid messages are being printed, it needs to print error message
//             //it is not an issue with the values for alpha, beta, zeta, or tripTime necessarily, so it will not print the cout about the type
//             notFirstSwitch = true; 
//         }
//         if (!skipPrint){ //if a message should be printed, goes through these statements
//             if (!notFirstSwitch){
//                 switch (i)
//                 {
//                 case ALPHA_OFFSET:
//                     cout << "alpha corresponded with the actual values for ";
//                     break;
//                 case BETA_OFFSET:
//                     cout << "beta corresponded with the actual values for ";
//                     break;
//                 case ZETA_OFFSET:
//                     cout << "zeta corresponded with the actual values for ";
//                     break;
//                 case TRIPTIME_OFFSET:
//                     cout << "trip time corresponded with the actual values for ";
//                     break;
//                 default:
//                     cout << "issue with verification function... ";
//                     noErrors = false; //something weird is going on with i which is also an error
//                     break;
//                 }
//             }
//             switch (whichMethod){
//             case 1:
//                 cout << "crossOver_wholeRandom";
//                 break;
//             case 2:
//                 cout << "crossOver_average";
//                 break;
//             case 3:
//                 cout << "crossOver_bundleVars";
//                 break;
//             default:
//                 cout << "incorrectly specified cross over method";
//                 noErrors = false; //if the method was not correctly specified, this is an error
//                 break;
//             }
//             if (utcConstants->thruster_type == thruster<double>::NO_THRUST){
//                 cout << " with no thruster" <<endl;
//             }
//             else {
//                 cout << " with a thruster" << endl;
//             }
//         }
//         if (!noErrors){ //if there were errors, does not bother to loop again because likely the other parameters will be broken 
//             return false;
//         }
//         parValIndex += 2; //must increase this at the end of each time through the loop so it's comparing to the correct numbers
//     }
//     return true;
// }

// //just makes checkReasonability shorter - takes in an offset and accesses a child's parameter corresponding to this offset
// double getParamStuff(const int & correspondingOffset, const Child& aChild){
//     if (correspondingOffset == ALPHA_OFFSET){ 
//         return aChild.startParams.alpha;
//     }
//     else if (correspondingOffset == BETA_OFFSET){
//         return aChild.startParams.beta;
//     }
//     else if (correspondingOffset == ZETA_OFFSET){
//         return aChild.startParams.zeta;
//     }
//     else if (correspondingOffset == TRIPTIME_OFFSET){
//         return aChild.startParams.tripTime;
//     }
//     else { //if there is an error with what is passed in, returns an error message and a weird, small negative
//         cout << "ERROR: undefined offset" << endl;
//         return -1e-06; //there was no specific reason for this number and it can be changed
//     }
// }

// void twentyAdultsPosAndSpeedDiffMade(bool printThings, std::vector<Adult>& allAdults, cudaConstants* utcConstants){
//     utcConstants->num_individuals = 20;

//     Child* genZero = new Child[utcConstants->num_individuals]; 
//     allAdults.clear();

//     //made 20 children with semi-random speed and position differences and tripTimes
//     //the tripTime has no actual correlation with the posDiff and speedDiff
//     //some of these tripTimes values are not correctly
//     genZero[0] = Child(40000000.0, 0.02, 0.0043); //about 1.27 years  
//     genZero[1] = Child(35000000.0, 0.048, 0.00234); //about 1.12 years 
//     genZero[2] = Child(41000000.0, 0.299, 0.0034); //about 1.30 years
//     genZero[3] = Child(32000000.0, 0.025, 0.00354); //about 1.01 years 
//     genZero[4] = Child(47000000.0, 0.122, 0.0034); //about 1.49 years 
//     genZero[5] = Child(43000000.0, 0.02, 0.0034); //about 1.36 years 
//     genZero[6] = Child(37400000.0, 0.012, 0.0042); //about 1.18 years 
//     genZero[7] = Child(38000000.0, 0.14, 0.00192); //about 1.20 years 
//     genZero[8] = Child(44300000.0, 0.053, 0.00414); //about 1.40 years 
//     genZero[9] = Child(45200000.0, 0.098, 0.00432); //about 1.43 years 
//     genZero[10] = Child(40000010.0, 0.02, 0.0043); //about 1.27 years  
//     genZero[11] = Child(35000020.0, 0.048, 0.00234); //about 1.12 years 
//     genZero[12] = Child(41000030.0, 0.299, 0.0034); //about 1.30 years
//     genZero[13] = Child(32000040.0, 0.025, 0.00354); //about 1.01 years 
//     genZero[14] = Child(47000050.0, 0.122, 0.0034); //about 1.49 years 
//     genZero[15] = Child(43000060.0, 0.02, 0.0034); //about 1.36 years 
//     genZero[16] = Child(37400070.0, 0.017, 0.0042); //about 1.18 years 
//     genZero[17] = Child(38000080.0, 0.14, 0.00192); //about 1.20 years 
//     genZero[18] = Child(44300090.0, 0.02, 0.0034); //about 1.40 years 
//     genZero[19] = Child(45200100.0, 0.299, 0.0034); //about 1.43 years

//     //all the children are converted intp Adults
//     for (int i = 0; i < utcConstants->num_individuals; i++){
//         allAdults.push_back(Adult(genZero[i]));
//     }

//     //ranks and gives distances to all the adults and then sorts them
//     giveRank(allAdults, utcConstants);
//     giveDistance(allAdults, utcConstants);
//     std::sort(allAdults.begin(), allAdults.end(), rankDistanceSort);

//     //if printThings is true, it will print all the individuals in rankDistanceSort order
//     if (printThings){
//         for (int i = 0; i < utcConstants->num_individuals; i++){
//             cout << "allAdults[" << i << "]: " << allAdults[i].unitTestingRankDistanceStatusPrint() << "- posDiff: " << allAdults[i].posDiff << ", speedDiff: " << allAdults[i].speedDiff << " & tripTime: " << allAdults[i].startParams.tripTime << endl;
//         }
//         cout << endl;
//     }

//     delete[] genZero; //cleans up dynamic memory
// }

// //checks to make sure children are being generated the way we think they are
// bool verifyChildrenFromCrossover(std::mt19937_64& rng, bool printThings, cudaConstants* utcConstants){
//     //holds the potential parents
//     std::vector<Adult> oldAdults;
//     //holds these potential parents based on whether or not they are duplicates or unique/original values
//     std::vector<Adult> parents;
//     std::vector<Adult> duplicates;
//     //this vector is necessary because now findDuplicates takes in two vectors, so since we only want to find the duplicates...
//     //  ...in oldAdults, we pass this empty vector in along with oldAdults
//     std::vector<Adult> empty;
//     //holds the children to be generated
//     Child* children = new Child[utcConstants->num_individuals]; 
//     Child occupiesLastSpace;
//     //starts the number of survivors at the entire population of individuals
//     utcConstants->survivor_count = utcConstants->num_individuals;
//     //holds the status of the code - whether or not we ran into any errors -> this will be the output
//     bool noErrors = true;

//     //the power of two which is dividing num_individuals to get the survivpr count
//     int exponent = 0;

//     //loops 4 times the number of survivors is n/(2^exponent) -> survivor_count = n (20), n/2 (10), n/4 (5), n/8 (2)
//     while (exponent < 4){
//         //ALL UNIQUE INDIVIDUALS    
//         //if it going to be printing information, it tell the user the survivor size and population size 
//         if (printThings){
//             cout << "Survivor population: " << utcConstants->survivor_count << " out of " << utcConstants->num_individuals << " potential parents" << endl;
//             cout << "With a pool of entirely unique individuals..." << endl;
//         }
//         bool goodThisTime = true; //initially assumes the test will pass
//         bool printParents = false; //set the print status to false because we do not need to print the 20 adults again before they are modified

//         //fills oldAdults with the 20 adults who only have posDiff, speedDiff, and tripTime uniquely assigned to them
//         twentyAdultsPosAndSpeedDiffMade(printParents, oldAdults, utcConstants);
//         empty.clear();

//         //assigns a duplicate status to all duplicates in oldAdults
//         findDuplicates(oldAdults, empty, utcConstants, utcConstants->anneal_initial); 

//         //a vector of unique posDiffs to replace the posDiff of any duplicate with
//         std::vector<double> posDiffs;
//         posDiffs.push_back(0.23);
//         posDiffs.push_back(0.153);
//         posDiffs.push_back(0.284);
//         posDiffs.push_back(0.032);
//         posDiffs.push_back(0.835);
//         posDiffs.push_back(0.556);
//         posDiffs.push_back(0.679);
//         posDiffs.push_back(0.332);
//         posDiffs.push_back(0.026);

//         //a vector of unique speedDiffs to replace the speedDiff of any duplicate with
//         std::vector<double> speedDiffs;
//         speedDiffs.push_back(0.009221);
//         speedDiffs.push_back(0.003865);
//         speedDiffs.push_back(0.008965);
//         speedDiffs.push_back(0.007342);
//         speedDiffs.push_back(0.00856);
//         speedDiffs.push_back(0.001356);
//         speedDiffs.push_back(0.0014368);
//         speedDiffs.push_back(0.008365);
//         speedDiffs.push_back(0.009243);
        
//         //make all the Adults in oldAdults unique individuals
//         for (int i = 0; i < oldAdults.size(); i++){
//             //if an individual is a duplicate, replace both its posDiff and speedDiff with different values to make it unique
//             if (oldAdults[i].errorStatus == DUPLICATE){
//                 oldAdults[i].posDiff = posDiffs[posDiffs.size()-1];
//                 posDiffs.pop_back();
//                 oldAdults[i].speedDiff = speedDiffs[speedDiffs.size()-1];
//                 speedDiffs.pop_back();
//             }
//         }

//         //ensure empty is actually empty, then tries to locate duplicates in oldAdukts
//         //there should no longer be any duplicates left
//         empty.clear();
//         findDuplicates(oldAdults, empty, utcConstants, utcConstants->anneal_initial);

//         //ensures all the adults in oldAdults are unique
//         for (int i = 0; i < oldAdults.size(); i++){
//             //if there is a duplicate, it prints an error message and will print oldAdults so the user can examine the issue
//             if (oldAdults[i].errorStatus == DUPLICATE){
//                 cout << i << endl;
//                 if (speedDiffs.size() == 0 || posDiffs.size() == 0){
//                     cout << "There were not enough unique speedDiffs and posDiffs to fix all the duplicate values or findDuplicates isn't working as expected" << endl;
//                 }
//                 else {
//                     cout << "Other unknown error having to do with proper posDiffs and speedDiffs not being assigned" << endl;
//                 }
//                 printParents = true;
//             }
//         }
        
//         //prints all the identifying information (including the error status for each Adult in oldAdults)
//         if (printParents){
//             for (int i = 0; i < utcConstants->num_individuals; i++){
//                 cout << "oldAdults[" << i << "]: " << oldAdults[i].unitTestingRankDistanceStatusPrint() << "- posDiff: " << oldAdults[i].posDiff << ", speedDiff: " << oldAdults[i].speedDiff << " & tripTime: " << oldAdults[i].startParams.tripTime << endl;
//             }
//             cout << endl;
//             printParents = false; //it will not do this large print next time 
//         }

//         //no mutations should occur because anneal_initial is set to 0, as is the default_mutation_chance
//         utcConstants->anneal_initial = 0;
//         utcConstants->default_mutation_chance = 0;

//         //set the generation at 1 because this is the first generation and that should have little influence on the generation
//         //the number of children to be generated using this method is num_individuals because there should be no duplicates at this point
//         generateChildrenFromCrossover(parents, children, utcConstants->num_individuals, rng, utcConstants->anneal_initial, 1, utcConstants);
//         goodThisTime = cfcAnswersMatchExpectations(occupiesLastSpace, utcConstants->num_individuals, children, utcConstants, parents, printThings);
//         //if the function that is meant to check if the answers match our expectations returns false, there are errors
//         if (!goodThisTime){
//             //identifies that there were errors if printThings is true
//             if (printThings){
//                 cout << "\n...generateChildrenFromCrossover() does not work as expected" << endl;
//             }
//             noErrors = false;
//         }

//         //A MIX OF UNIQUE INDIVIDUALS AND DUPLICATES
//         //identifies the next test being done
//         if (printThings){
//             cout << "In a pool that is a mix of unique individuals and duplicates...";
//         }
//         //refills oldAdults with 11 unique individuals and 9 duplicates 
//         twentyAdultsPosAndSpeedDiffMade(printParents, oldAdults, utcConstants);
//         empty.clear();
//         findDuplicates(oldAdults, empty, utcConstants, utcConstants->anneal_initial);
        
//         //generates individuals based on solely the original/unique individuals from the population
//         //generation is set to two to distinguish individuals created in this second call from those created in the first
//         generateChildrenFromCrossover(oldAdults, children, utcConstants->num_individuals, rng, utcConstants->anneal_initial, 2, utcConstants);
//         goodThisTime = cfcAnswersMatchExpectations(occupiesLastSpace, utcConstants->num_individuals, children, utcConstants, oldAdults, printThings);
//         //if the function that is meant to check if the answers match our expectations returns false, there are errors
//         if (!goodThisTime){
//             //identifies that there were errors if printThings is true
//             if (printThings){
//                 cout << "\n...generateChildrenFromCrossover() does not work as expected" << endl;
//             }
//             noErrors = false;
//         }

//         //ALL DUPLICATES
//         //identifies the next test being done
//         if (printThings){
//             cout << "In a pool that is entirely duplicates of one individual...";
//         }
//         //fill oldAdults entirely with duplicates 
//         for (int i = 0; i < oldAdults.size(); i++){
//             oldAdults[i] = Adult(Child(3.8e+07, 0.014, 0.00192));
//             //prints all the identifying information (including the error status for each Adult in oldAdults) if printParents is true
//             if (printParents){
//                 cout << "oldAdults[" << i << "]: " << oldAdults[i].unitTestingRankDistanceStatusPrint() << "- posDiff: " << oldAdults[i].posDiff << ", speedDiff: " << oldAdults[i].speedDiff << " & tripTime: " << oldAdults[i].startParams.tripTime << endl;
//             }
//         }
//         //if it printed all these
//         if (printParents){
//             cout << endl;
//             printParents = false; //it will not do this large print next time 
//         }
//         //locates the duplicates
//         empty.clear();
//         findDuplicates(oldAdults, empty, utcConstants, utcConstants->anneal_initial);

//         //generates individuals based on solely the original/unique individual from oldAdults
//         //can the code handle only one parent? If so, what does it do?
//         //generation is set to three to distinguish individuals created in this third call from those created in the first and second
//         generateChildrenFromCrossover(oldAdults, children, utcConstants->num_individuals - duplicates.size(), rng, utcConstants->anneal_initial, 3, utcConstants);
//         goodThisTime = cfcAnswersMatchExpectations(occupiesLastSpace, utcConstants->num_individuals - duplicates.size(), children, utcConstants, oldAdults, printThings);
//         cout << "got past third generate and check" << endl;
//         //if the function that is meant to check if the answers match our expectations returns false, there are errors
//         if (!goodThisTime){
//             //identifies that there were errors if printThings is true
//             if (printThings){
//                 cout << "\n...generateChildrenFromCrossover() does not work as expected" << endl;
//             }
//             noErrors = false;
//         }

//         utcConstants->survivor_count /= 2; //Halves the survivor population
//         exponent++; //increases the exponent
//     }
//     delete[] children; //deallocates dynamic memory
    
//     return noErrors; //returns whether or not the tests were successful
// }

// // Compares the results from verifyChildrenFromCrossover with the correct values
// bool cfcAnswersMatchExpectations(const Child & endSpot, const int & numChildren, const Child* childrenGenerated, const cudaConstants* utcConstants, std::vector<Adult> & parents, bool printThings){
//     bool noErrors = true; //initially assumes that all the values passed in will meet our expectations

//     //if there numChildren does not fill the entire childrenGenerated array, ensures that it did not overwrite a value that was not supposed to change
//     if (numChildren < utcConstants->num_individuals){
//         if(childrenGenerated[numChildren].birthday != endSpot.birthday || childrenGenerated[numChildren].posDiff != endSpot.posDiff || childrenGenerated[numChildren].speedDiff != endSpot.speedDiff){
//             noErrors = false;
//             cout << "The last child was overwritten" << endl;
//         }
//     }
    
//     int currSet = 0; //the set of children this is examining
//     while(currSet*crossoverChildrenCount < numChildren){
//         bool correctParents = false; //true if we can be fairly sure a Child had two unique and valid parents
//         bool ableToCheckBothParents = false; //this is only false if we have one more child than would make a complete set (7, 13, 19, etc)
//         int prnt1 = -1, prnt2 = -1; //represents the index of a Child's parents
//         //if a pair of parents have generated 3 or more children
//         if (currSet*crossoverChildrenCount + fullAvgOffset < numChildren){ 
//             ableToCheckBothParents = true; //we will be able to identify two unique parents
//             for (int i = 0; i < parents.size(); i++){
//                 //a parent must have a partner to generate a Child, so it is counted as an error if a parent generates a Child by itself
//                 for (int j = i+1; j < parents.size(); j++){ 
//                     if ((parents[i].startParams.tripTime + parents[j].startParams.tripTime)/2 == childrenGenerated[currSet*crossoverChildrenCount + fullAvgOffset].startParams.tripTime){
//                         correctParents = true;
//                         prnt1 = i;
//                         prnt2 = j;
//                         //prints the tripTimes of the parents and of any children they've made
//                         if (printThings){
//                             cout << std::fixed << std::setprecision(1) << "Parent 1's tripTime: " << parents[prnt1].startParams.tripTime << "\t\t Parent 2's tripTime: " << parents[prnt2].startParams.tripTime << endl;

//                             int kids = 0; //how many children this parent set had

//                             //loops until it has printed the tripTime for every Child a couple made
//                             while(currSet*crossoverChildrenCount + kids < numChildren && kids < crossoverChildrenCount){
//                                 cout << std::fixed << std::setprecision(1) << "Child " << kids+1 << "'s tripTime: " << childrenGenerated[currSet*crossoverChildrenCount+kids].startParams.tripTime;
//                                 //prints two Children's tripTimes per line, so only ends the line every other time
//                                 if (kids%2 == 1){
//                                     cout << endl;
//                                 }   
//                                 else {
//                                     cout << "\t\t ";
//                                     if (currSet*crossoverChildrenCount + kids+1 == numChildren){
//                                         cout << endl;
//                                     }
//                                 }
//                                 kids++; //increments kids and checks if the parents created another one
//                             }  
//                             cout << endl; //prints an extra blank line between every output
//                         }
//                     }
//                 }
//             }
//         }
//         //if a pair of parents generated only 1 or 2 children
//         else{
//             //attempts to find a parent whose tripTime values match the first child's tripTimes values because this child is generated by the crossOver_wholeRandom method
//             for (int i = 0; i < parents.size(); i++){
//                 //if a parent's tripTime matches the first Child's tripTime, it was a parent of this incomplete set of children
//                 if (parents[i].startParams.tripTime == childrenGenerated[currSet*crossoverChildrenCount].startParams.tripTime){
//                    prnt1 = i;
//                 }
//                 //if there is a second child that can be checked, we may be able to identify the other parent
//                 if (currSet*crossoverChildrenCount + 1 < numChildren){
//                     ableToCheckBothParents = true; //there is a chance we can identify both parents
//                     //if a parent's tripTime matches this second child's tripTime it is also assumed to be a parent of the incomplete set
//                     if (parents[i].startParams.tripTime == childrenGenerated[currSet*crossoverChildrenCount + 1].startParams.tripTime){
//                         prnt2 = i;
//                     }
//                 }
//             }
//             //if printThings is true, it will print the parents' and their children's tripTimes
//             if (printThings){
//                 cout << std::fixed << std::setprecision(1) << "Parent 1's tripTime: " << parents[prnt1].startParams.tripTime << "\t\t Parent 2's tripTime: ";
//                 //if the code identified two distinct parnets, print each of their tripTimes and their Children's
//                 //to be able to identify two unique parents this way, there must be two Children so prints both
//                 if (ableToCheckBothParents && prnt1 != prnt2){
//                     cout << std::fixed << std::setprecision(1) << parents[prnt2].startParams.tripTime << endl;
//                     cout << std::fixed << std::setprecision(1) << "Child 1's tripTime: " << childrenGenerated[currSet*crossoverChildrenCount].startParams.tripTime;
//                     cout << std::fixed << std::setprecision(1) << "\t\t Child 2's tripTime: " << childrenGenerated[currSet*crossoverChildrenCount +1].startParams.tripTime << endl;
//                 }
//                 //if it somehow fails to detect to distinct parents, but has two different Children (this should never happen)
//                 //prints both Children's tripTimes and identifies that the second parent's tripTime is unknown
//                 else if (ableToCheckBothParents){
//                     cout << "UNKNOWN" << endl;
//                     cout << std::fixed << std::setprecision(1) << "Child 1's tripTime: " << childrenGenerated[currSet*crossoverChildrenCount].startParams.tripTime;
//                     cout << std::fixed << std::setprecision(1) << "\t\t Child 2's tripTime: " << childrenGenerated[currSet*crossoverChildrenCount +1].startParams.tripTime << endl;
//                 }
//                 //if there was only one child, it could not identify the second parent
//                 //prints a message that it could not identify the second parent's tripTime 
//                 //and prints the known Child's tripTime
//                 else {
//                     cout << "UNKNOWN" << endl;
//                     cout << std::fixed << std::setprecision(1) << "Child 1's tripTime: " << childrenGenerated[currSet*crossoverChildrenCount].startParams.tripTime << endl;
//                 }
//                 cout << endl; //prints an extra blank line at the end
//             }
//         }
//         //if prnt1 was given a value
//         if (prnt1 > -1){
//             //it will identify whether it only checked one or both parents
//             if (ableToCheckBothParents){
//                 //if parent 2 has been identified and is unique, it has found both parents
//                 if (prnt2 > -1 && prnt2 != prnt1){
//                     correctParents = true;
//                 }
//                 //otherwise, there was an error
//                 else{
//                     correctParents = false;
//                 }
//             }
//             //if only one parent was identified, it cannot conclusively say whether or not the set would have the correct parents...
//             //  ...since tripTimes was the only variable that was set
//             else{
//                 cout << "Unable to determine if children were correctly generated from parents" << endl;
//                 correctParents = true; 
//             }
//         }
//         //if prnt1 is not even identified, it has not correctly identified the parents
//         else {
//             correctParents = false;
//         }
//         //if it did not correctly identify the parents, then there are errors
//         if (!correctParents){
//             if (printThings){
//                 cout << "Could not correctly identify the parent" << endl;
//             }
//             noErrors = false;
//         }
//         currSet++; //goes on to checking the next set of children
//     }
//     return noErrors; //if all children had their parents identified and they seem good, this returns true
// }

// bool testingCopyOfNewGen(bool printThings, std::mt19937_64 & rng, cudaConstants* utcConstants){
//     //originally assumes that all the functions are working the way we expect them to
//     bool allWorking = true;
//     //starts the number of survivors and individuals at 20 so we can use all the duplicates to make children
//     utcConstants->num_individuals = 20;
//     utcConstants->survivor_count = 10;
//     //ensures anneal_initial is set as well as the default and duplicate mutation chances
//     utcConstants->anneal_initial = 0.01;
//     utcConstants->default_mutation_chance = 0.0;

//     //holds the potential parents
//     std::vector<Adult> oldAdults;
//     std::vector<Adult> newAdults;

//     //holds the newAdults once they have been converted back into children for verification
//     Child* newChildren = new Child[utcConstants->num_individuals]; 

//     //sets the generation to 1 because generation should not have a major impact on any of these function
//     int generation = 1;

//     //the hard-coded false is just supposed to stop this from printing all 20 Adults again
//     twentyAdultsPosAndSpeedDiffMade(false, oldAdults, utcConstants);

//     //UTCopyOfNewGen returns the number of individuals that should have been generated by crossOver
//     int endSpotHere = UTCopyOfNewGeneration(oldAdults, newAdults, utcConstants->anneal_initial, generation, rng, utcConstants);

//     //turns the Adult that is the first that should have been generated by mutation into endSpot
//     //this Child is just created because cfcAnswersMatchExpectations needs it
//     Child endSpot(newAdults[endSpotHere].startParams, newAdults[endSpotHere].posDiff, newAdults[endSpotHere].speedDiff, newAdults[endSpotHere].progress, newAdults[endSpotHere].avgParentProgress, newAdults[endSpotHere].birthday);

//     //the two verification functions (cfcAnswersMatchExpectations and properMutations need an array of Children as their inputs)
//     convertBackToChildren(newAdults, newChildren, utcConstants);

//     //if the children from crossover were not correctly generted, it makes the values that will be returned false 
//     if (!cfcAnswersMatchExpectations(endSpot, endSpotHere, newChildren, utcConstants, oldAdults, printThings)){
//         allWorking = false;
//     }

//     //if the test worked as expected, returns TRUE
//     return allWorking;
// }

// //converts Adults to Children so they can be put into the existing verification algorithms
// void convertBackToChildren(std::vector<Adult>& newAdult, Child* newChildren, const cudaConstants* utcConstants){
//     //both newAdults and newChildren should be num_individuals size
//     for (int i = 0; i < utcConstants->num_individuals; i++){
//         //uses a unit testing Child constructor to create Children with pretty much all the same elements as the Adult has
//         newChildren[i] = Child(newAdult[i].startParams, newAdult[i].posDiff, newAdult[i].speedDiff, newAdult[i].progress, newAdult[i].avgParentProgress, newAdult[i].birthday);
//         //if the posDiffs are not the same that is a good sign that something really isn't working with the Child constructor
//         if (newChildren[i].posDiff != newAdult[i].posDiff || newChildren[i].speedDiff != newAdult[i].speedDiff){
//             cout << "This needs a better constructor " << endl;
//         }
//     }
// }

// //copied from genetic_algorithm.cpp on 6/29/22 at 5:59PM and then removed all the callRK and callRK dependent things
// //made it return an int representing where the children from the unique individuals ended and the duplicates started
// int UTCopyOfNewGeneration(std::vector<Adult> & oldAdults, std::vector<Adult> & newAdults, const double & annealing, const int & generation, std::mt19937_64 & rng, const cudaConstants* utcConstants){
//     //Vector that will hold the adults who are potential parents
//     //The criteria for being a parent is being in the top survivor_count number of adults in the oldAdult pool
//     std::vector<Adult> parents; 

//     parents.clear();

//     //separates oldAdults into parents (full of best N/4 individuals) 
//     fillParents(oldAdults, parents, generation, utcConstants);

//     //Import number of new individuals the GPU needs to fill its threads
//     //Create a newChildren function to fill in with children generated from oldAdults
//     Child* newChildren = new Child[utcConstants->num_individuals]; 

//     //uses parents to make children for the next generation
//     makeChildren(parents, newChildren, annealing, generation, rng, utcConstants);
    
//     //Now that the children have been simulated, convert the children into adults
//     //First, it will calculate the right pos and speed diffs for the children
//     //This will also put the converted children into the newAdults vector
//     convertToAdults(newAdults, newChildren, utcConstants); 

//     //Free the pointer's memory
//     delete[] newChildren; 
// }

// //unit testing creating an entire generation of individuals
// bool firstFullGen(std::mt19937_64& rng, cudaConstants * utcConstants, bool printThings){
//     utcConstants->num_individuals = 10;
//     utcConstants->survivor_count = 3;

//     bool noErrors = true;

//     Child* genZero = new Child[utcConstants->num_individuals]; 
//     elements<double> elems; //makes a set of default elements
//     coefficients<double> coeffs = coefficients<double>(); //makes a set of default coefficients

//     rkParameters<double>* paramsForIndividuals = new rkParameters<double>[utcConstants->num_individuals];
    
//     //creates parameters for making child with unique tripTimes, but default elements and default coefficients 
//     paramsForIndividuals[0] = rkParameters<double>(40000000.0, elems, coeffs); //about 1.27 years  
//     paramsForIndividuals[1] = rkParameters<double>(35000000.0, elems, coeffs); //about 1.12 years 
//     paramsForIndividuals[2] = rkParameters<double>(41000000.0, elems, coeffs); //about 1.30 years
//     paramsForIndividuals[3] = rkParameters<double>(32000000.0, elems, coeffs); //about 1.01 years 
//     paramsForIndividuals[4] = rkParameters<double>(47000000.0, elems, coeffs); //about 1.49 years 
//     paramsForIndividuals[5] = rkParameters<double>(43000000.0, elems, coeffs); //about 1.36 years 
//     paramsForIndividuals[6] = rkParameters<double>(37400000.0, elems, coeffs); //about 1.18 years 
//     paramsForIndividuals[7] = rkParameters<double>(38000000.0, elems, coeffs); //about 1.20 years 
//     paramsForIndividuals[8] = rkParameters<double>(44300000.0, elems, coeffs); //about 1.40 years 
//     paramsForIndividuals[9] = rkParameters<double>(45200000.0, elems, coeffs); //about 1.43 years 
    
//     //turns the rkParameters into children to make up genZero
//     for (int i = 0; i < utcConstants->num_individuals; i++){
//         //all these children were created in the first generation and have no parents, so the parent 
//         genZero[i] = Child(paramsForIndividuals[i], utcConstants, 1, 0);
//     }

//     std::vector<Adult> parents;
//     firstGeneration(genZero, parents, utcConstants); //genZero and is turned into parents using firstGeneration which has been unit tested previously

//     //the parents are sorted so they can be selected to be parents for a new generation
//     giveRank(parents, utcConstants); 
//     std::sort(parents.begin(), parents.end(), rankSort); 
//     giveDistance(parents, utcConstants); 
//     std::sort(parents.begin(), parents.end(), rankDistanceSort);

//     std::vector<Adult> youngGen;

//     //creates young a young generation with an annealing rate of 0.0 and the generation is set to 1 (the random 0.0 and 1)
//     newGeneration(parents, youngGen, 0.0, 1, rng, utcConstants);

//     //verifies the children generated are as expected
//     noErrors = verifyFullGen(youngGen, parents, utcConstants, printThings);

//     //if the two vectors do not have the same size, there is definitely an issue
//     if (youngGen.size() != parents.size() || parents.size() != utcConstants->num_individuals){
//         noErrors = false;
//     }

//     //deallocates dynamic memory
//     delete[] genZero;
//     delete[] paramsForIndividuals;

//     //if the first test passed successfully, tries making a lot a children from a small number of survivors
//     if (noErrors){
//         cout << "Test successful - trying creating a large number number of children from a small number of adults " << endl;
//         noErrors = makeManyChildren(rng, youngGen, parents, utcConstants, printThings);
//     }

//     //returns noErrors
//     return noErrors;
// }

// // A function that is used to verify that firstFullGen is working correctly
// //Made this earlier before we started treating duplicates differently and so long as duplicates are handled the same as regular numbers
// //     this test should still work
// bool verifyFullGen(std::vector<Adult>& youngGen, std::vector<Adult>& possParents, const cudaConstants * utcConstants, bool printThings){

//     cout << utcConstants->survivor_count << "; " << utcConstants->num_individuals << endl;
//     bool noErrors = true;

//     std::vector<int> parentIndexSets; //vector that holds the indices of sets of parents

//     //int setsOfChildren = (utcConstants->num_individuals)/childrenPerSet;
//     int setNum = -1; //starts at negative 1 because incremented on the first run

//     int indexVecSize = 0;
//     double tripTimeTol = 10e-6;

//     //it checks which parents are being used and prints things if printThings is set to true
//     for (int i = 0; i < utcConstants->num_individuals; i++){
//         //when it has gone through a whole set generated by one set of parents, it increase setNum
//         if (i%crossoverChildrenCount == 0){ 
//             setNum++;
//         }
//         printThings = true;
//         if (printThings){ //prints all the parents and children at each index
//             if (i < possParents.size()){
//                 if (i < utcConstants->survivor_count){ //puts *| |* around any parent that has been selected to have children
//                     cout << "*|";
//                 }
//                 cout << "parents[" << i << "]'s tripTime: ";
//                 cout << possParents[i].startParams.tripTime;
//                 if (i < utcConstants->survivor_count){
//                     cout << "|*";
//                 }
//                 else {
//                     cout << "\t";
//                 }
//             }
//             else if (i > utcConstants->survivor_count){
//                 cout << "                               ";
//             }   
//             cout << "\t\t";

//             cout << "youngGen[" << i << "]'s tripTime: ";
//             cout << youngGen[i].startParams.tripTime << endl;
//         }
//         //the averaging function is the best way to determine which two parents created a child
//         //none of the survivor count individuals should have the same average
//         if (i == fullAvgOffset + crossoverChildrenCount*setNum){
//             indexVecSize += 2; //each time this occurs, two parents should be added to the parentIndexSets vector
//             //if the average tripTime is equal to the copied value of parent, just one parent created the child and it is a clone
//             if (youngGen[i].startParams.tripTime == youngGen[6*setNum].startParams.tripTime){ 
//                 cout << "The parent did not have a partner and probably created clones" << endl;
//                 noErrors = false; //this is an error if a parent clones itself
//                 for (int j = 0; j < utcConstants->survivor_count; j++){
//                     if (youngGen[i].startParams.tripTime == possParents[j].startParams.tripTime){
//                         cout << "Set " << setNum << " created by parent[" << j << "] and parent[" << j << "] " << endl;
//                         parentIndexSets.push_back(j);
//                         parentIndexSets.push_back(j);
//                     }
//                 }
//             }
//             else{
//                 //otherwise, it searches different pairs of potentially chosen parents to determine which parents produced this offspring
//                 for (int j = 0; j < utcConstants->survivor_count; j++){
//                     for (int k = j+1; k < utcConstants->survivor_count; k++){
//                         //checks the two parents
//                         if (youngGen[i].startParams.tripTime <= ((possParents[j].startParams.tripTime + possParents[k].startParams.tripTime)/2.0 + tripTimeTol) && youngGen[i].startParams.tripTime >= ((possParents[j].startParams.tripTime + possParents[k].startParams.tripTime)/2.0 - tripTimeTol)){
//                             cout << "Set " << setNum << " created by parent[" << j << "] and parent[" << k << "] " << endl;
//                             parentIndexSets.push_back(j);
//                             parentIndexSets.push_back(k);
//                         }
//                     }
//                 }
//             }
//             //if the two parents have not been added, this is an issue and it prints error messages, then pushes back -1s so we will be able to detect if this happens again
//             if(indexVecSize != parentIndexSets.size()){
//                 noErrors = false;
//                 cout << indexVecSize << " vs " << parentIndexSets.size() << endl;
//                 cout << "Unknown parents for set #" << setNum << " of children" << endl;
//                 parentIndexSets.push_back(-1);
//                 parentIndexSets.push_back(-1);

//             }
//         }
//     }
//     return noErrors;
// }

// //tries making children from just a couple of survivors - not currently working as far as I know...
// bool makeManyChildren(std::mt19937_64& rng, std::vector<Adult>& youngGen, std::vector<Adult>& possParents, cudaConstants * utcConstants, bool printThings){
//     bool noErrors = true; 

//     utcConstants->num_individuals = 65;
//     utcConstants->survivor_count = 7;

//     newGeneration(possParents, youngGen, 0.0, 1, rng, utcConstants);
//     //verifies 
//     noErrors = verifyFullGen(youngGen, possParents, utcConstants, false);

//     if (!noErrors){
//         cout << "Problem with 65 individuals generated from 7 parents" << endl;
//     }
//     else if (printThings){
//         cout << "Generated 65 individuals from 7 parents" << endl;
//     }

//     utcConstants->survivor_count = 5;
//     newGeneration(possParents, youngGen, 0.0, 1, rng, utcConstants);
//     //verifies 
//     noErrors = verifyFullGen(youngGen, possParents, utcConstants, false);

//     if (!noErrors){
//         cout << "Problem with 65 individuals generated from 5 parents" << endl;
//     }
//     else if (printThings){
//         cout << "Generated 65 individuals from 5 parents" << endl;
//     }

//     return noErrors;
// }

// bool checkUTMutateMask(){
//     //set up the mask to be mutated
//     bool * mutateMask= new bool[OPTIM_VARS];
//     //rng based on time seed
//     std::mt19937_64 rng(time(NULL));
//     //mutation rate from config
//     double mutation_rate = 0.75;
//     //call the UT fuction
//     UTmutateMask(rng, mutateMask, mutation_rate);

//     //print the results
//     cout << "Mask: " << endl;
//     for(int i = 0; i < OPTIM_VARS; i++){
//         cout << i << " M: " << mutateMask[i] << endl;
//     }

//     return true;
    
// }

// void UTmutateMask(std::mt19937_64 & rng, bool * mutateMask, double mutation_rate){
    
//     for (int i = 0; i < OPTIM_VARS; i++) {
//         //Reset mask
//         mutateMask[i] = false;
//     }
    
//     // counter to make sure in the unlikely event that all genes are being mutated the code doesn't
//     // get stuck in infinite loop looking for a false spot in the array
//     int geneCount = 0; //How many gene's have we changed
//     // Set a gene to mutate if a randomized values is less than mutation_rate, repeating everytime this is true
//     while ((static_cast<double>(rng()) / rng.max()) < mutation_rate && geneCount < OPTIM_VARS) {
//         bool geneSet = false; // boolean to flag if a gene was successfully selected to be set as true
//         int index; // index value that will be assigned randomly to select a gene
//         while (geneSet != true) {
//             index = rng() % OPTIM_VARS;
//             // If the randomly selected gene hasn't been already set to mutate, set it and flag the geneSet to true to break out of the loop
//             if (mutateMask[index] == false) {
//                 mutateMask[index] = true;
//                 geneSet = true;
//                 geneCount++;
//             }
//         }

//     }
//     //make sure this matched how many are set to true
//     cout << "Final geneCount: " << geneCount << endl;
// }
