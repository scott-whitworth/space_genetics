#include <math.h>
#define _USE_MATH_DEFINES // for use of M_PI

// Default constructor
Adult::Adult(){
    rank = INT_MAX; //might change later?
    distance = -1;
    rarity = INT_MAX;
    associatedRefPoint = 0;
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Compare two individuals by their positional difference values, used in standard sort to have array contain lowest posDiff individual at start
// input: two individuals
// output: returns true if personB has a higher positional difference than personA or true/false based on the adults' error statuses
bool LowerPosDiff(Adult& personA, Adult& personB) {
    if(personA.errorStatus !=  VALID && personA.errorStatus != DUPLICATE){
        return false;
    }
    else if(personB.errorStatus !=  VALID && personB.errorStatus != DUPLICATE){
        return true;
    }
    else if (personA.posDiff < personB.posDiff) {
        return true;
    }
    else {
        return false;
    }

}

// Compare two individuals by their velocity difference values, used in standard sort to have array contain highest speedDiff individual at start
// input: two individuals
// output: returns true if personB has a lower velocity difference than personA
// bool HigherSpeedDiff(const Adult& personA, const Adult& personB) {
//     if(personA.errorStatus !=  VALID && personA.errorStatus != DUPLICATE){
//         return false;
//     }
//     else if(personB.errorStatus !=  VALID && personB.errorStatus != DUPLICATE){
//         return true;
//     }
//     else if (personA.speedDiff > personB.speedDiff) {
//         return true;
//     }
//     else {
//         return false;
//     }
// }

// Compare two individuals by their velocity difference values, used in standard sort to have array contain lowest speedDiff individual at start
// input: two individuals
// output: returns true if personA has a lower velocity difference than personB
bool LowerSpeedDiff(const Adult& personA, const Adult& personB) {
    if(personA.errorStatus !=  VALID && personA.errorStatus != DUPLICATE){
        return false;
    }
    else if(personB.errorStatus !=  VALID && personB.errorStatus != DUPLICATE){
        return true;
    }
    else if (personA.speedDiff < personB.speedDiff) {
        return true;
    }
    else {
        return false;
    }
}

// Compare two adults by their orbit positional difference values
// input: two adults
// output: returns true if personB has a lower orbit pos difference than personA
bool LowerOrbitPosDiff(Adult& personA, Adult& personB) {
    if(personA.errorStatus !=  VALID && personA.errorStatus != DUPLICATE){
        return false;
    }
    else if(personB.errorStatus !=  VALID && personB.errorStatus != DUPLICATE){
        return true;
    }
    else if (personA.orbitPosDiff < personB.orbitPosDiff) {
        return true;
    }
    else {
        return false;
    }
}

// Compare two adults by their orbit speed difference values
// input: two adults
// output: returns true if personB has a lower orbit speed difference than personA
bool LowerOrbitSpeedDiff(Adult& personA, Adult& personB) {
    if(personA.errorStatus !=  VALID && personA.errorStatus != DUPLICATE){
        return false;
    }
    else if(personB.errorStatus !=  VALID && personB.errorStatus != DUPLICATE){
        return true;
    }
    else if (personA.orbitSpeedDiff < personB.orbitSpeedDiff) {
        return true;
    }
    else {
        return false;
    }    
}

// Compare two individuals by their spent fuel values, used in standard sort to have array contain lowest fuel spent individual at start
// input: two individuals
// output: returns true if personA has a smaller amount of fuel used than personB
bool LowerFuelSpent(const Adult& personA, const Adult& personB) {
    if(personA.errorStatus !=  VALID && personA.errorStatus != DUPLICATE){
        return false;
    }
    else if(personB.errorStatus !=  VALID && personB.errorStatus != DUPLICATE){
        return true;
    }
    else if (personA.fuelSpent < personB.fuelSpent) {
        return true;
    }
    else {
        return false;
    }
}

// Compare two individuals by their triptime values, used in standard sort to have array contain lowest triptime individual at start
// input: two individuals
// output: returns true if personA has a lower triptime than personB
bool LowerTripTime(const Adult& personA, const Adult& personB) {
    if(personA.errorStatus !=  VALID && personA.errorStatus != DUPLICATE){
        return false;
    }
    else if(personB.errorStatus !=  VALID && personB.errorStatus != DUPLICATE){
        return true;
    }
    else if (personA.startParams.tripTime < personB.startParams.tripTime) {
        return true;
    }
    else {
        return false;
    }
}

// Compare two individuals by their minMarsDist values, used in standard sort to have array contain lowest minMarsDist individual at start
// input: two individuals
// output: returns true if personA has a lower minMarsDist than personB
bool LowerMarsDist(const Adult& personA, const Adult& personB) {
    if(personA.errorStatus !=  VALID && personA.errorStatus != DUPLICATE){
        return false;
    }
    else if(personB.errorStatus !=  VALID && personB.errorStatus != DUPLICATE){
        return true;
    }
    else if (personA.minMarsDist < personB.minMarsDist) {
        return true;
    }
    else {
        return false;
    }
}

// Compare two individuals by their orbithChange values, used in standard sort to have array contain highest orbithChange individual at start
// input: two individuals
// output: returns true if personA has a lower orbithChange than personB
bool LowerOrbitHChange(const Adult& personA, const Adult& personB){
        if(personA.errorStatus !=  VALID && personA.errorStatus != DUPLICATE){
        return false;
    }
    else if(personB.errorStatus !=  VALID && personB.errorStatus != DUPLICATE){
        return true;
    }
    else if (personA.orbithChange < personB.orbithChange) {
        return true;
    }
    else {
        return false;
    }
}

//compares two individuals and sorts based on which individual is closer to 1
//WARNING: this may be implemented in output.cpp in printFinalGen() when it may not need to be (use for test reasons)
//input: two individuals
//output: returns true if personA is closer to 1 than personB
bool bestProgress(const Adult& personA, const Adult& personB) {
    if(personA.errorStatus !=  VALID && personA.errorStatus != DUPLICATE){
        return false;
    }
    if(personB.errorStatus !=  VALID && personB.errorStatus != DUPLICATE){
        return true;
    }
    if (personA.progress > personB.progress) {
        return true;
    }
    else {
        return false;
    }

}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Compare two individuals to see if the first individual dominates the second individual
//Returns true if personA dominates personB.
//returns false if personA does not dominate personB.
bool dominationCheck(Adult& personA, Adult& personB, const cudaConstants* cConstants) {
    //Tracks if person A is better than person B in at least some compared parameter
    //  Domination requires only one parameter to be better for A to be better than B (as long as the others are at least equal)
    //  Thus, bool will start as false, so any parameter where A is better than B makes the bool true
    bool aParamBetter = false;

    //For loop will run through the objectives to check all necessary parameters against eachother
    //  If it is found that A has worse parameters than B (when tolerance has not been met), it will return false since A needs to be at least equal to dominate
    for (int i = 0; i < cConstants->missionObjectives.size(); i++) {

        //Create temp doubles to store the adults' parameters that will be compared
        double aParam = personA.getParameters(cConstants->missionObjectives[i]);
        double bParam = personB.getParameters(cConstants->missionObjectives[i]);
        
        //First check determines if the goal is to favor higher values or if it is to favor lower values
        //  This is done by looking at the value of the objective's goal
        //  If it is below 0, it means that the program seeks to minimize
        //  Otherwise, if it is above 0, the program seeks to maximize
        //  A value of 0 means there is an error with the goal
        //This is needed because the signs of the comparisons need to be changed
        // if (cConstants->missionObjectives[i].goal < 0) {
            //While not converged ...
            if ( !(aParam < cConstants->missionObjectives[i].dominationThreshold
                && bParam < cConstants->missionObjectives[i].dominationThreshold) ) {
                
                //Means they haven't met the parameters
                //First check to see if A's parameter is less than B's parameter (within a tolerance)
                if (aParam < (bParam - cConstants->missionObjectives[i].equateTolerance)) {
                    //Since A is better in at least one parameter, set the aParamBetter flag to true
                    aParamBetter = true;
                }
                //If here, it means that the adults haven't met the threshold and that A's parameter is not better than B's
                //Now there is a check for if A's parameters are worse/larger than B's (within a tolerance)
                else if (aParam > (bParam + cConstants->missionObjectives[i].equateTolerance)) {
                    //Since a parameter for A is worse than B, A doesn't dominate B, so return false
                    return false; 
                }
                //If here (without triggering the first if statement), it must mean that A's and B's params are the same (within a tolerance)
                //Now, the program will continue on to the next objective (if there are any)
            }
            // else: both parameters are past the dominationThreshold
            //       this means: we are ignoring that parameter for domination check
            
        // }
        //The maximization objectives operate very similarly to minimization, except that the signs are flipped
        // else if (cConstants->missionObjectives[i].goal > 0) {
        //     //Check to see if see if the parameters are both above/better than the tolerence, which means that no comparison will be necessary
        //     if ( !(aParam > cConstants->missionObjectives[i].dominationThreshold
        //         && bParam > cConstants->missionObjectives[i].dominationThreshold) ) {
                
        //         //Means they haven't met the parameter's tolerance
        //         //First check to see if A's parameter is greater than B's parameter (within a tolerance)
        //         if (aParam > (bParam + cConstants->missionObjectives[i].equateTolerance)) {
        //             //Since A is better in at least one parameter, set the aParamBetter flag to true
        //             aParamBetter = true;
        //         }
        //         //If here, it means that the adults haven't met the threshold and that A's parameter is not better than B's
        //         //Now there is a check for if A's parameters are worse/smaller than B's (within a tolerance)
        //         else if (aParam < (bParam - cConstants->missionObjectives[i].equateTolerance)) {
        //             //Since a parameter for A is worse than B, A doesn't dominate B, so return false
        //             return false; 
        //         }
        //         //If here (without triggering the first if statement), it must mean that A's and B's params are the same (within a tolerance)
        //         //Now, the program will continue on to the next objective (if there are any)
        //     }
        // }
        // else {
        //     std::cout << "\n_-_-_-_-_-_-_-_-_-Error Identifying Parameter Goal_-_-_-_-_-_-_-_-_-\n";
        // }
    }
    //The program have compared the parameters for the adults
    //  Now, if the program is here it means that A is at least as good as B (for all parameters)
    //  But, this could mean that either A truly dominates B, or they could be duplicates
    //  So, it is necessary to return the status of aParamBetter
    //      If it's true, it means one of A's parameters are better than B, and since the other parameters are at least as good, A dominates B
    //      If it is fale, all of A's values are the same as B, which means they are duplicates
    return aParamBetter; 
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Compare two individuals by their rank
//WARNING: Using this function to sort individuals will sort them by rank, but within those ranks they will be sorted by the last method used to sort them.
//         For example, if all the individuals are rank 1, sorting them using this method will do nothing. 
//input: two individuals
//output: if person A's rank is lower than person B's rank, return true
bool rankSort(const Adult& personA, const Adult& personB) {
    if(personA.errorStatus != VALID && personA.errorStatus != DUPLICATE){ //if personA has nan values or other errors they are set as worse than other adults (even other adults with errors)
        return false;
    }
    else if (personB.errorStatus != VALID && personB.errorStatus != DUPLICATE){ //if personA is not a nan, but personB is, personA is better
        return true;
    }
    if (personA.rank < personB.rank) {
        return true;
    }
    else {
        return false;
    }
}

//Compare two individuals by their rank and distance
//input: two individuals
//output: if person A's rank is lower than person B's rank, return true
//        if person A and person B have the same rank and person A has a greater distance than person B, return true
//Sorts the whole pool from lowest to highest rank. Individuals of the same rank are sorted from highest to lowest distance
bool rankDistanceSort(const Adult& personA, const Adult& personB) {
    if (personA.errorStatus != VALID && personA.errorStatus != DUPLICATE){
        return false;
    }
    else if (personB.errorStatus != VALID && personB.errorStatus != DUPLICATE){ //if personA is not a nan, but personB is, personA is better
        return true;
    }
    if(personA.rank < personB.rank){
        return true;
    }
    else if (personA.rank == personB.rank && personA.distance > personB.distance){
        return true;
    }
    else {
        return false;
    }  

}

//Compare two individuals by their rank and rarity
//input: two individuals
//output: if person A's rank is lower than person B's rank, return true
//        if person A and person B have the same rank and person A has a better (lower) rarity than person B, return true
//Sorts the whole pool from lowest to highest rank. Individuals of the same rank are sorted from lowest (best) to highest (worst) rarity score
bool rankRaritySort(const Adult& personA, const Adult& personB) {
    if (personA.errorStatus != VALID && personA.errorStatus != DUPLICATE){
        return false;
    }
    else if (personB.errorStatus != VALID && personB.errorStatus != DUPLICATE){ //if personA is not a nan, but personB is, personA is better
        return true;
    }
    if(personA.rank < personB.rank){
        return true;
    }
    else if (personA.rank == personB.rank && personA.rarity < personB.rarity){
        return true;
    }
    else {
        return false;
    }  

}

bool duplicateCheck(const Adult& personA, const Adult& personB, const cudaConstants* cConstants){

    //For loop will iterate through all objectives/parameters
    for (int i = 0; i < cConstants->missionObjectives.size(); i++) {
        //First check to see if A's parameter for this objective is lower than B's within a tolerance
        if (personA.getParameters(cConstants->missionObjectives[i]) < (personB.getParameters(cConstants->missionObjectives[i]) - cConstants->missionObjectives[i].equateTolerance)) {
            //Means that one of A's parameters doesn't equal B, so return false
            return false;
        }
        //Now there is a need to check if A's parameter is larger than B's
        else if (personA.getParameters(cConstants->missionObjectives[i]) > (personB.getParameters(cConstants->missionObjectives[i]) + cConstants->missionObjectives[i].equateTolerance)) {
            //If here, it means that A's parameter is larger than B's, so they aren't duplicates 
            return false;
        }
    }

    //If the program is here, it means that it has gone through all of A's and B's parameters without returning false
    //  This measn that it didn't find any differences in the parameters, so they are duplicates
    //  Thus, we need to return true
    return true; 

}

//Find the duplicates within the imported vector
void findDuplicates (std::vector<Adult>& newAdults, std::vector<Adult>& oldAdults, const cudaConstants* cConstants, const double& currentAnneal) {
    //reset the status of all the duplicates
    //This makes sure that previous duplicates get another look if their copies have been changed
    std::vector<Adult> tempAllAdults;
    int newAdultsSize = newAdults.size();
    int oldAdultsSize = oldAdults.size();

    
    for(int i = 0; i < oldAdults.size(); i++){
        if(oldAdults[i].errorStatus == DUPLICATE){
            oldAdults[i].errorStatus = VALID;
        }
        tempAllAdults.push_back(oldAdults[i]);
    }
    for(int i = 0; i < newAdults.size(); i++){
        if(newAdults[i].errorStatus == DUPLICATE){
            newAdults[i].errorStatus = VALID;
        }
        tempAllAdults.push_back(newAdults[i]);
    }
    oldAdults.clear();
    newAdults.clear();
    //loop through all the adults
    for (int i = 0; i < tempAllAdults.size(); i++){
        //i+1 so it doesn't check itself or past indexes
        for(int j = i+1; j < tempAllAdults.size(); j++){
            //only true if it is both a duplicate and has not been previous marked as a duplicate
            // [j].duplicate check is for the second time an Adult is flagged as a duplicate
            //Checks for a valid error status, so duplicates can be reset to valid without worry of overriding other error statuses
            if(duplicateCheck(tempAllAdults[i], tempAllAdults[j], cConstants) && tempAllAdults[j].errorStatus == VALID) {
                tempAllAdults[j].errorStatus = DUPLICATE;
            }
        }
    }

    for(int i = 0; i < oldAdultsSize; i++){
        oldAdults.push_back(tempAllAdults[i]);
    }

    for(int i = oldAdultsSize; i < (newAdultsSize + oldAdultsSize); i++){
        newAdults.push_back(tempAllAdults[i]);
    }
}
#ifdef UNITTEST //this should be defined in unit testing

bool duplicateCheckTest(const Adult& personA, const Adult& personB){

    //true if A posDiff "equals" B posDiff
    bool APosEqualsB = false;
    //true if A speedDiff "equals" B speedDiff
    bool ASpeedEqualsB = false;

    //tolerances used to determine the range of values considered equal
    //these are both currently set to 1e-14 AU, I don't think these need to be modified 
    double posTolerance = 1e-14;
    double speedTolerance = 1e-14;

    //True is A posdiff is equal to B posDiff +- posTolerance
    if ((personA.posDiff < personB.posDiff + posTolerance) && (personA.posDiff > personB.posDiff - posTolerance)){
        APosEqualsB = true;
    }
    //True is A speeddiff is equal to B speedDiff +- speedTolerance
    if ((personA.speedDiff < personB.speedDiff + speedTolerance) && (personA.speedDiff > personB.speedDiff - speedTolerance)){
        ASpeedEqualsB = true;
    }

    //if they are both true then we found a duplicate
    if(APosEqualsB == true && ASpeedEqualsB == true){
        return true;
    }else{//a dulplicate was not found
        
        return false;
    }

    return false;
}

bool dominationCheckTest(Adult& personA, Adult& personB, int missionType){
    
    // //Is true if A is at least equally as good as B for all objectives
    // bool AisEqual = false;
    // //Is true if A is better than B for at least one objective
    // bool AisBetter = false;
    // //tolerances used to determine the range of values considered equal
    // //these are both currently set to 1e-14 AU, I don't think these need to be modified 
    // //this tolerance is about 0.0015m, and I don't think we can go lower?

    // double posTolerance = 5;
    // double speedTolerance = 5;
    // //double posTolerance = 0.0;
    // //double speedTolerance = 0.0;
    // //true if A posDiff "equals" B posDiff
    // bool APosEqualsB = false;
    // //true if A speedDiff "equals" B speedDiff
    // bool ASpeedEqualsB = false;

    // //std::cout << "\n\nPERSON A TEST: posDiff = " << personA.posDiff << ", speedDiff = " << personA.speedDiff; 
    // //std::cout << "\n\nPERSON B TEST: posDiff = " << personB.posDiff << ", speedDiff = " << personB.speedDiff; 

    // //True is A posdiff is equal to B posDiff +- posTolerance
    // if ((personA.posDiff < personB.posDiff + posTolerance) && (personA.posDiff > personB.posDiff - posTolerance)){
    //     APosEqualsB = true;
    // }
    // //True is A speeddiff is equal to B speedDiff +- speedTolerance
    // if ((personA.speedDiff < personB.speedDiff + speedTolerance) && (personA.speedDiff > personB.speedDiff - speedTolerance)){
    //     ASpeedEqualsB = true;
    // }
    // //If the mission type is a rendezvous, A's speed diff needs to be lower for it to be equal or better
    // if (missionType == Impact) {
    //     /*//Check if A's posDiff is approximately/better than B's posDiff and check the same for speed
    //     //If so, A is equal to B
    //     if ((personA.posDiff < personB.posDiff || APosEqualsB) && (personA.speedDiff > personB.speedDiff || ASpeedEqualsB)) {
    //         AisEqual = true; 
    //     }  
    //     //If A has a better posDiff or speedDiff than B, then A is better than B
    //     if (personA.posDiff < personB.posDiff || personA.speedDiff > personB.speedDiff) {
    //         AisBetter = true;
    //     }*/
    //     //For impact, only optimizing for one opbjective
    //     if ((personA.posDiff < personB.posDiff || APosEqualsB)) {
    //         AisEqual = true; 
    //     }  
    //     //If A has a better posDiff or speedDiff than B, then A is better than B
    //     if (personA.posDiff < personB.posDiff) {
    //         AisBetter = true;
    //     }
    // }
    // //If the mission type is a impact, A's speed diff needs to be higher than B's for it to be equal or better
    // else {
    //     //If A.posDiff is approximately/better than B.posDiff, and A.speedDiff is approximately/better than B.speedDiff, then A is equal to B.
    //     if ((personA.posDiff < personB.posDiff || APosEqualsB) && (personA.speedDiff < personB.speedDiff || ASpeedEqualsB)) {
    //         AisEqual = true;
    //     }
    //     //If A has a better posDiff or speedDiff than B, then A is better than B
    //     if (personA.posDiff < personB.posDiff || personA.speedDiff < personB.speedDiff){
    //         AisBetter = true;
    //     }
        
    // }

    // //A Dominates B
    // //A only dominates B if:
    //     //A.posDiff - posTolerance < B.posDiff <= A.posDiff, and A.speedDiff < B.speedDiff
    //     //A.posDiff < B.posDiff, and A.speedDiff - speedTolerance < B.speedDiff <= A.speedDiff
    //     //A.posDiff < B.posDiff and A.speedDiff < B.speedDiff (better in every way)
    // if (AisEqual && AisBetter){
    //     //std::cout << "\n\nDomination check returns TRUE";
    //     return true;
    // }
    // //A and B are codominant if:
    //     //A.posDiff < B.posDiff & A.speedDiff > B.speedDiff
    //     //A.posDiff > B.posDiff & A.speedDiff < B.speedDiff
    //     //A.posDiff = B.posDiff & A.speedDiff = B.speedDiff (Unlikely since they're doubles, unless A & B are clones or genetically similar)
    // //B dominates A if any of the conditions for "A dominates B", when reversed, are true.
    // else {
    //     //std::cout << "\n\nDomination check returns FALSE";
     return false;
    // }
}

const std::string Adult::unitTestingRankDistanceStatusPrint(){
    std::string allInfo = "(" + std::to_string(rank) + "," + std::to_string(distance);
    if (errorStatus == VALID){
        allInfo += ", VALID)";
    }
    else {
        allInfo += ", Invalid - error type)";
    }
    return allInfo;
}

int Adult::getRank(){
    return rank;
}

double Adult::getDistance(){
    return distance;
}

#endif //unit testing endif