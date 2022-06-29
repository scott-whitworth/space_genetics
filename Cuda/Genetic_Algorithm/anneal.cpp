//Function that will adjust the annneal based on the previous anneal and if there was a change in the best individual
void changeAnneal (const std::vector<Adult>& oldAdults, const cudaConstants* cConstants, double & currentAnneal, int & oldestBirthday, double & dRate){
    
    //Grab the current best (based on rank-distance) cost of this generation
    double curCost = oldAdults[0].cost;
    double previousAnneal = currentAnneal;
    //Caluclate the new current anneal
    //It will be a linear decrease from the initial/max anneal based on how well the current cost is compared to the cost threshold

    //How we calculate the anneal depends on the mission type
    //  Our current formula is anneal_max * (1 - (Mission goals / cost))
    //  The section in the parentheses should equal 0 when the cost has been met
    //  Thus, the number of mission goals are different between the missions, so we need to check what type of mission it is to calculate the anneal
    if (cConstants -> missionType == Impact) 
    {
        currentAnneal = cConstants->anneal_initial * (1 - (curCost));
    }
    else 
    {
        if(oldestBirthday > 500 && oldestBirthday % 10 == 0){
            dRate = dRate/1.0;
        }

        currentAnneal = (cConstants->anneal_initial * (1 - pow((curCost), 0.1)))*dRate;
        if(currentAnneal > previousAnneal){
            currentAnneal = previousAnneal;
        }
    }
    

    //Check to make sure that the current anneal does not fall below the designated minimum amount
    if (currentAnneal < cConstants->anneal_final)
    {
        currentAnneal = cConstants->anneal_final;
    }
}

//TODO: Do we need this anymore now that cost is calulated for each individual?
//Function that will calculate this generation's best adult's cost
/*
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
*/