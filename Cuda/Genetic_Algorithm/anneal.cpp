//TODO: Figure out what to replace posDiff (what used to be cost)
//Function that will adjust the annneal based on the previous anneal and if there was a change in the best individual
void changeAnneal (const std::vector<Adult>& oldAdults, const cudaConstants* cConstants, double & currentAnneal, int & oldestBirthday,  double & previousBestPosDiff, int & generation, const double & posTolerance, double & dRate){
    
    //Calculate the current cost for this generation
    double curCost = calculateCost(oldAdults, cConstants);
    
    //Caluclate the new current anneal
    //It will be a linear decrease from the initial/max anneal based on how well the current cost is compared to the cost threshold

    //How we calculate the anneal depends on the mission type
    //  Our current formula is anneal_max * (1 - (Mission goals / cost))
    //  The section in the parentheses should equal 0 when the cost has been met
    //  Thus, the number of mission goals are different between the missions, so we need to check what type of mission it is to calculate the anneal
    if (cConstants -> missionType == Impact) 
    {
        currentAnneal = cConstants->anneal_initial * (1 - (Impact / curCost));
    }
    else 
    {
        if(oldestBirthday > 500 && oldestBirthday % 10 == 0){
            dRate = dRate/1.3;
        }
        currentAnneal = (cConstants->anneal_initial * (1 - pow((Rendezvous / curCost), 0.18)))*dRate;
    }
    

    //Check to make sure that the current anneal does not fall below the designated minimum amount
    if (currentAnneal < cConstants->anneal_final)
    {
        currentAnneal = cConstants->anneal_final;
    }
    

    /*
    // Scaling anneal based on proximity to tolerance
    // Far away: larger anneal scale, close: smaller anneal
    if (cConstants->missionType == Impact) {
        //Impact is only based on posDiff, so proximity-based annealing only relies on how close posDiff is to tolerance.
        new_anneal = currentAnneal * (1 - (posTolerance / oldAdults[0].posDiff));
    }

    else if (cConstants->missionType == Rendezvous) {
        if (posTolerance < oldAdults[0].posDiff){ 
            //TODO: decide what we want to do with this annealing   
            //Exponentially changing annealing TODO: this annealing and the impact annealing are hardly getting adjusted
            new_anneal = currentAnneal * (1 - pow(posTolerance / oldAdults[0].posDiff,2.0));
            if (new_anneal < cConstants->anneal_final){
                new_anneal = cConstants->anneal_final; //Set a true minimum for annealing
            }
        }
    }
    
    /*
    //Process to see if anneal needs to be adjusted
    // If generations are stale, anneal drops
    Adult currentBest;
    // Compare current best individual to that from CHANGE_CHECK (50) many generations ago.
    // If they are the same, change size of mutations
    if (static_cast<int>(generation) % (cConstants->change_check) == 0) { 
        currentBest = oldAdults[0];
        // checks for anneal to change
        // previousBest starts at 0 to ensure changeInBest = true on generation 0
        if ( !(changeInBest(previousBestPosDiff, currentBest, dRate)) ) { 
            //this ensures that changeInBest never compares two zeros, thus keeping dRate in relevance as the posDiff lowers
            if (trunc(currentBest.posDiff/dRate) == 0) { //posDiff here used to be cost
                while (trunc(currentBest.posDiff/dRate) == 0) { //posDiff here used to be cost
                    dRate = dRate/10; 
                }
                std::cout << "\nnew dRate: " << dRate << std::endl;
            }
            // If no change in BestIndividual across generations, reduce currentAnneal by anneal_factor while staying above anneal_min
            //reduce anneal_min
            anneal_min = cConstants->anneal_initial*exp(-sqrt(posTolerance/oldAdults[0].posDiff)*generation);
            if (anneal_min < cConstants->anneal_final){
                anneal_min = cConstants->anneal_final;//Set a true minimum for annealing
            }

            //Rendezvous mission uses anneal_min, Impact does not
            if(cConstants->missionType == Impact) {
                currentAnneal = currentAnneal * cConstants->anneal_factor;
            }
            else if (cConstants->missionType == Rendezvous){
                currentAnneal = (currentAnneal * cConstants->anneal_factor > anneal_min)? (currentAnneal * cConstants->anneal_factor):(anneal_min);
            }
            std::cout << "\nnew anneal: " << currentAnneal << std::endl;              
        }

        previousBestPosDiff = currentBest.posDiff; //posDiff here used to be cost
    }
    */
}
//TODO: What should we do with old methods now that the new one works well?
//tests if the best value has changed much in the since the last time this was called
bool changeInBest(double previousBestPosDiff, const Adult & currentBest, double distinguishRate) {
    //truncate is used here to compare doubles via the distinguguishRate, to ensure that there has been relatively no change.
        if (trunc(previousBestPosDiff/distinguishRate) != trunc(currentBest.posDiff/distinguishRate)) {
            return true;
        }
        else { 
            return false;
        }
}

//Function that will calculate this generation's best adult's cost
double calculateCost(const std::vector<Adult> & oldAdults, const cudaConstants* cConstants){
    //Create the cost double that will be returned
    double cost; 

    //Variables will be used to calculate cost
    //      They both will specifically be used to calculate the adult's diff / goal diff
    //      This prevents a excellent speed diff causing the adult's diff / goal diff to be less than 1
    //          If this was the case, there is the potential that the cost goal could be met even if the position diff hadn't hit it's goal
    double rel_posDiff, rel_speedDiff;

    //How to calculate the cost will depend on the mission type
    //In general, cost will be (for each mission parameter, the difference/the goal) - the number of mission parameters
    //This will mean a cost of 0 or below signifies that the individual has hit the mission goals

    //For impacts, the cost only depends on posDiff
    //  NOTE: While the RD sorting favors high speed for impacts, convergence ultimately comes down to posDiff
    if (cConstants -> missionType == Impact) {
        //The goal is position threshold and the current status is the best individual's posDiff

        //Check to see if the rel_posDiff would be less than the threshold
        if (oldAdults[0].posDiff > cConstants->pos_threshold) {
            //If not, calculate how close the diff is to the goal
            rel_posDiff = oldAdults[0].posDiff / cConstants->pos_threshold; 
        }
        else {
            //If so, set the rel_posDiff to 1, signifying that the goal has been met
            rel_posDiff = 1.0;
        }

        //The coast is the rel_posDiff minus the number of mission goals (1 in this case)
        cost = rel_posDiff;
    }
    //For rendezvous, the cost depends on both posDiff and speedDiff
    else {
        //Similarly to impact, calculate how far the best adult's speed and position diffs are away from the goal and subtract by the number of mission goals
        
        //First, same as above either set rel_posDiff to how close the best adult's posDiff is to the goal or to 1, depending on if the posDiff is beating the goal
        if (oldAdults[0].posDiff > cConstants->pos_threshold) {
            rel_posDiff = oldAdults[0].posDiff / cConstants->pos_threshold;
        }
        else {
            rel_posDiff = 1.0;
        }
        //Repeat the same process in calculating rel_posDiff for rel_speedDiff
        if (oldAdults[0].speedDiff > cConstants->speed_threshold) {
            rel_speedDiff = oldAdults[0].speedDiff / cConstants->speed_threshold;
        }
        else { 
            rel_speedDiff = 1.0;
        }

        //The cost the addition of each difference between the goals and pod/speed diffs minus the rendezvous mission goal count (2)
        //This means a flight that meets all goals will have a cost of 0, since the rel pos/speed diffs would be set to 1
        cost = (rel_posDiff + rel_speedDiff);
    }
    
    //Return the calculated cost
    return cost; 
}