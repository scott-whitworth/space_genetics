//returns true if the first generation is generated
bool firstParentsTest(){
    std::vector<Adult> parents;
    Child* theChildren = new Child[genSize];
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

    //TODO: I assme the following VVV is a 'key' to the above, but it might be reasonable to have a rationale for the layout of the above ^^^

    //Ranks based on my net dominations (times it dominates others - times it is dominated by others) -> dominations based on output of dominationCheckTest
    //Rank 1: (20, 70) [NET DOMINATIONS: 6]
    //Rank 2: (150, 20); (110, 40); (30, 90) [NET DOMINATIONS: 3]
    //Rank 3: (180, 30); (20, 120) [NET DOMINATIONS: 1]
    //Rank 4: (120, 90) [NET DOMINATIONS: -1]
    //Rank 5: (50, 120) [NET DOMINATIONS: -2]
    //Rank 6: (340, 90) [NET DOMINATIONS: -6]
    //Rank 7: (220, 970) [NET DOMINATIONS-7]

    UTfirstGeneration(theChildren, parents);
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

//Simpified unit test version of first generation -> does not callRK
void UTfirstGeneration(Child* initialChildren, std::vector<Adult>& oldAdults){
    UTconvertToAdults(oldAdults, initialChildren);

}

//Simplified unit test version of convert to adults
void UTconvertToAdults(std::vector<Adult> & newAdults, Child* newChildren){
    //Iterate through the newChildren vector to add them to the newAdult vector
    //iterates until num_individuals as that is the size of newChildren
    for (int i = 0; i < genSize; i++)
    {
        //Fill in the newAdults vector with a new adult generated with a completed child
        newAdults.push_back(Adult(newChildren[i]));
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
    int k = 0; //a number k that represents tthe rank that will be assigned to an individual (k is a bad name for it, that was a holdover from when the below was a for loop not a while loop)
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

bool firstFullGen(){
    
}

void UTnewGen(std::vector<Adult> & oldAdults, std::vector<Adult> & newAdults, const double & annealing, const int & generation, std::mt19937_64 & rng){
    //Import number of new individuals the GPU needs to fill its threads
    //Create a newChildren function to fill in with children generated from oldAdults
    Child* newChildren = new Child[genSize]; 

    //Create a mask to determine which of the child's parameters will be inherited from which parents
    std::vector<int> mask;

    for (int i = 0; i < OPTIM_VARS; i++)
    {
        mask.push_back(AVG); //TODO: Why are we setting this to AVG?
    }

    //Array that determines which parents will be paired up to generate children
    //it will contain indexes of adults that are deemed to be parents
    //it will be shuffled to generate random parent pairings
    //NOTE: this will only be effective after oldAdults has been sorted (with the best Adults at the front)
    //int parentPool [(cConstants->survivor_count)] = {0};
    std::vector<int> parentPool;

    //Fill the parentPool index with the number of indexes desired for survivors
    for (int i = 0; i < cConstants->survivor_count; i++)
    {
        parentPool.push_back(i);
    }

    //Shuffle the parentPool mapping array
    //This will generate a random list of indicies, which can then be used to map to random Adults within oldAdult's survivor range to pair parents
    //This will effectively generate random parent pairs
    std::shuffle(parentPool.begin(), parentPool.end(), rng); 

    //Int that tracks the number of new individuals that have been generated
    //Used primarily to make sure that no children are overwritten in newChildren; indexing in generatechildrenPair, should equal N when finished
    //Starts as 0 beacuse no new children have been generated
    int numNewChildren = 0;

    //Int that tracks the index of parents; determines which parents are passed in to the generate children function
    //Will also help track if the oldAdults array needs to be reshuffled to generate new parent pairs
    int parentIndex = 0;

    //While loop will make sure to generate children based on how many GPU threads are available
    //This will result in children being generated until there are enough to fill the gpu
    //Note: at the time of writing, the number of threads is tied to num_individuals
    while (numNewChildren < cConstants->num_individuals) {
        //TODO: You do have a potential here that is mapped in generateChildrenPair
        //      If you generate 3 pairs (so 6 children), and your num_individuals is not evenly divisible by 6, then you are going to have a memory issue
        //      There is also a potential problem is this algorithm ends up needing more than int parentNum = cConstants->num_individuals/4; pairs
        //      I like the start, but think deeply about how each of these elements are going to be moving around

        //Generate a pair of children based on the random cross over mask method from a pair of parents
        //Generate the base mask
        crossOver_wholeRandom(mask, rng);

        //TODO: Clean up these comments. If you ever copy/paste anything you are probably doing it wrong

        //Generate a pair of children based on the mask
        generateChildrenPair(oldAdults[parentPool[parentIndex]], oldAdults[parentPool[parentIndex+1]], newChildren, mask, annealing, rng, numNewChildren, generation, cConstants);

        //Generate a pair of children from a pair of parents based on the average mask method
        //Generate the base mask
        crossOver_average(mask);

        //Generate a pair of children based on the mask
        generateChildrenPair(oldAdults[parentPool[parentIndex]], oldAdults[parentPool[parentIndex+1]], newChildren, mask, annealing, rng, numNewChildren, generation, cConstants);

        //Generate a pair of children from a pair of parents based on the bundled random mask method
        //Generate the base mask
        crossOver_bundleVars(mask, rng);

        //Generate a pair of children based on the mask
        generateChildrenPair(oldAdults[parentPool[parentIndex]], oldAdults[parentPool[parentIndex+1]], newChildren, mask, annealing, rng, numNewChildren, generation, cConstants);

        //TODO: My main interest is in getting rid of this second bundleVars crossover. If you feel you are ready, I would remove it
        //For a second time, generate a pair of children feom a pair of parents based on the bundled random mask method
        //Generate the base mask
        crossOver_bundleVars(mask, rng);

        //Generate a pair of children based on the mask
        generateChildrenPair(oldAdults[parentPool[parentIndex]], oldAdults[parentPool[parentIndex+1]], newChildren, mask, annealing, rng, numNewChildren, generation, cConstants);

        //Iterate through the shuffled section of the oldAdults array
        //Add two to parentIndex to account for the variable tracking pairs of parents, not just one parent's index
        parentIndex += 2;

        //Check to see if parents from outside the selected survivor/shuffled range
        //This will make sure that only the best adults have the potential to become parents
        if (parentIndex >= (cConstants->survivor_count)) {
            //Reset the parent index to start from the beginning of the shuffled section of oldAdults
            parentIndex = 0;

            //Re-shuffle the best/survivor section of oldAdults to make sure there are new pairs of parents
            //This will ideally not harm genetic diversity too much, even with the same set of parents
            std::shuffle(parentPool.begin(), parentPool.end(), rng);

            //TODO: Kind of weird, but I like it.
            //      Ok, so the rationale is that you only want to pull from the 'best' individuals
            //      You don't want to pull from anywhere in oldAdults, because the bottom 3/4 of oldAdults is deemed to be 'bad'
            //      Thus you will only ever consider the 'top' 1/4 of oldAdults, enforced by parentNum

            //TODO: It sounds like we need to make parentNum a configured value. I could see us playing around with modifying the size of the parent pool
            //      1/4 is kind of arbitrary, it might be interesting to see what happens if we kept 1/8 or some fixed number / percent
        }
    }

    std::cout << "\n\n_-_-_-_-_-_-_-_-TEST: POST CHILDREN CREATION_-_-_-_-_-_-_-_-\n\n";

    double timeInitial = 0;
    double calcPerS = 0;

    // Each child represents an individual set of starting parameters 
        // GPU based runge kutta process determines final position and velocity based on parameters
    //Will fill in the final variables (final position & speed, posDiff, speedDiff) for each child
    //TODO: get rid of magic numbers, first 0 is timeInitial, second 0 is calcPerS, the are both set as 0 before being passed in to callRK within optimize.cu
        //Perhaps we should just import cCOnstants and newChildren into callRk, since most of the arguments come from cConstants regardless
    callRK(cConstants->num_individuals, cConstants->thread_block_size, newChildren, timeInitial, (cConstants->orbitalPeriod / cConstants->GuessMaxPossibleSteps), cConstants->rk_tol, calcPerS, cConstants); 

    //Determine the status and diffs of the simulated children
    setStatusAndDiffs(newChildren, cConstants); 

    std::cout << "\n\n_-_-_-_-_-_-_-_-TEST: POST RK_-_-_-_-_-_-_-_-\n\n";

    //Now that the children have been simulated, convert the children into adults
    //This will also put the converted children into the newAdults vector
    convertToAdults(newAdults, newChildren, cConstants); 

    std::cout << "\n\n_-_-_-_-_-_-_-_-TEST: POST CONVERSION_-_-_-_-_-_-_-_-\n\n";

    //Free the pointer's memory
    delete[] newChildren; 
}