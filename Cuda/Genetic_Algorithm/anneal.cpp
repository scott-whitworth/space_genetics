//Function that will adjust the annneal based on the previous anneal and if there was a change in the best individual
void changeAnneal (const std::vector<Adult>& oldAdults, const cudaConstants* cConstants, double & currentAnneal, int & oldestBirthday, double & dRate){
    
    //Grab the current best (based on rank-distance) progress of this generation
    double curProgress = oldAdults[0].progress;
    double previousAnneal = currentAnneal;

    //Caluclate the new current anneal
    //It will be adjusted based on a power scale

    //How we calculate the anneal depends on the mission type
    //  Our current formula is anneal_initial * (1 - curProgress^0.1) * dRate
    //  dRate is used to help move anneal along when there isn't a ton of progress being made
    //  The section in the parentheses should equal 1 when the goals has been met
    //  However, to make sure the anneal doesn't go too low, we make sure that it doen't go below the min anneal
    //  Thus, the number of mission goals are different between the missions, so we need to check what type of mission it is to calculate the anneal

    //TODO: discuss if we need the distinction between the missions here since there is a distinction when calculating progress
    if (cConstants -> missionType == Impact) 
    {
        currentAnneal = cConstants->anneal_initial * (1 - (curProgress));
    }
    else 
    {
        if(oldestBirthday > 500 && oldestBirthday % 10 == 0){
            dRate /= 1.0;
        }

        currentAnneal = (cConstants->anneal_initial * (1 - pow((curProgress), 0.1)))*dRate;
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