#include <math.h>
#define _USE_MATH_DEFINES // for use of M_PI

// Default constructor
Adult::Adult(){
    rank = INT_MAX; //might change later?
    distance = -1;
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Compare two individuals by their positional difference values, used in standard sort to have array contain lowest posDiff individual at start
// input: two individuals
// output: returns true if personB has a higher positional difference than personA
bool LowerPosDiff(Adult& personA, Adult& personB) {
    if(personA.errorStatus !=  VALID && personA.errorStatus != DUPLICATE){
        return false;
    }
    if(personB.errorStatus !=  VALID && personB.errorStatus != DUPLICATE){
        return true;
    }
    if (personA.posDiff < personB.posDiff) {
        return true;
    }
    else {
        return false;
    }

}

// Compare two individuals by their velocity difference values, used in standard sort to have array contain highest speedDiff individual at start
// input: two individuals
// output: returns true if personB has a lower velocity difference than personA
bool HigherSpeedDiff(const Adult& personA, const Adult& personB) {
    if (personA.speedDiff > personB.speedDiff) {
        return true;
    }
    else {
        return false;
    }
}

// Compare two individuals by their velocity difference values, used in standard sort to have array contain lowest speedDiff individual at start
// input: two individuals
// output: returns true if personA has a lower velocity difference than personB
bool LowerSpeedDiff(const Adult& personA, const Adult& personB) {
    if (personA.speedDiff < personB.speedDiff) {
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
    
    //Is true if A is at least equally as good as B for all objectives
    bool AisEqual = false;
    //Is true if A is better than B for at least one objective
    bool AisBetter = false;
    //TODO: Might want to consider modifying tolerances 
    //tolerances used to determine the range of values considered equal
    //these are both currently set to 1e-14 AU, I don't think these need to be modified 
    //this tolerance is about 0.0015m, and I don't think we can go lower?

    //TODO:Might want to consider deleting most of this function

    //pos and speed variables used for comparing adutls to eachother
    double posTolerance = cConstants->posDominationTolerance;
    double speedTolerance = cConstants->speedDominationTolerance;
    double pos_threshold = cConstants-> posDominationThreshold;
    double speed_threshold = cConstants-> speedDominationThreshold;

    //true if A posDiff "equals" B posDiff
    bool APosEqualsB = false;
    //true if A speedDiff "equals" B speedDiff
    bool ASpeedEqualsB = false;

    //True is A posdiff is equal to B posDiff +- posTolerance
    if (((personA.posDiff < personB.posDiff + posTolerance) && (personA.posDiff > personB.posDiff - posTolerance)) || (personA.posDiff < pos_threshold && personB.posDiff < pos_threshold)){
        APosEqualsB = true;
    }
    //True is A speeddiff is equal to B speedDiff +- speedTolerance
    if (((personA.speedDiff < personB.speedDiff + speedTolerance) && (personA.speedDiff > personB.speedDiff - speedTolerance)) || (personA.speedDiff < speed_threshold && personB.speedDiff < speed_threshold)){
        ASpeedEqualsB = true;
    }
    //If the mission type is a rendezvous, A's speed diff needs to be lower for it to be equal or better
    if (cConstants -> missionType == Impact) {
        //For impact, only optimizing for one objective
        if ((personA.posDiff < personB.posDiff || APosEqualsB)) {
            AisEqual = true; 
        }  
        //If A has a better posDiff or speedDiff than B, then A is better than B
        if (personA.posDiff < personB.posDiff) {
            AisBetter = true;
        }
    }
    //If the mission type is a impact, A's speed diff needs to be higher than B's for it to be equal or better
    else {
        //If A.posDiff is approximately/better than B.posDiff, and A.speedDiff is approximately/better than B.speedDiff, then A is equal to B.
        if ((personA.posDiff < personB.posDiff || APosEqualsB) && (personA.speedDiff < personB.speedDiff || ASpeedEqualsB)) {
            AisEqual = true;
        }
        //If A has a better posDiff or speedDiff than B, then A is better than B
        if (personA.posDiff < personB.posDiff || personA.speedDiff < personB.speedDiff){
            AisBetter = true;
        }
        
    }

    //A Dominates B
    //A only dominates B if:
        //A.posDiff - posTolerance < B.posDiff <= A.posDiff, and A.speedDiff < B.speedDiff
        //A.posDiff < B.posDiff, and A.speedDiff - speedTolerance < B.speedDiff <= A.speedDiff
        //A.posDiff < B.posDiff and A.speedDiff < B.speedDiff (better in every way)
    if (AisEqual && AisBetter){
        //std::cout << "\n\nDomination check returns TRUE";
        return true;
    }
    //A and B are codominant if:
        //A.posDiff < B.posDiff & A.speedDiff > B.speedDiff
        //A.posDiff > B.posDiff & A.speedDiff < B.speedDiff
        //A.posDiff = B.posDiff & A.speedDiff = B.speedDiff (Unlikely since they're doubles, unless A & B are clones or genetically similar)
    //B dominates A if any of the conditions for "A dominates B", when reversed, are true.
    else {
        //std::cout << "\n\nDomination check returns FALSE";
        return false;
    }
    /*
    */
}

//!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Compare two individuals by their rank
//WARNING: Using this function to sort individuals will sort them by rank, but within those ranks they will be sorted by the last method used to sort them.
//         For example, if all the individuals are rank 1, sorting them using this method will do nothing. 
//input: two individuals
//output: if person A's rank is lower than person B's rank, return true
bool rankSort(const Adult& personA, const Adult& personB){
    //TODO: is this a good system to check validity first?
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

//Compare two adults based on their distance, sorting the lowest distances first
//  This function will be used to detect duplicates within mutateAdults in ga_crossover
//input:  PersonA - First adult to be compared
//        PersonB - Second adult to be compared
//output: True if personA's distance is less than personB's or if personB's status isn't valid
//        Fale if personB's distance is less than personA's or if personA's status isn't valid
//              Note: personA's status is checked before personB's, so if neither person is valid, it will return false
//  NOTE: This function assumes that the distance for both adults have already been calculated
bool lowerDistanceSort(const Adult& personA, const Adult& personB) {
    if(personA.errorStatus != VALID && personA.errorStatus != DUPLICATE){ //if personA has nan values or other errors they are set as worse than other adults (even other adults with errors)
        return false;
    }
    else if (personB.errorStatus != VALID && personB.errorStatus != DUPLICATE){ //if personA is not a nan, but personB is, personA is better
        return true;
    }
    if (personA.distance < personB.distance) {
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

bool duplicateCheck(const Adult& personA, const Adult& personB, const cudaConstants* cConstants, const double& currentAnneal){

    //true if A posDiff "equals" B posDiff
    bool APosEqualsB = false;
    //true if A speedDiff "equals" B speedDiff
    bool ASpeedEqualsB = false;

    //tolerances used to determine the range of values considered equal
    //these are both currently set to 1e-14 AU, I don't think these need to be modified 
    double posTolerance = cConstants->pos_threshold*currentAnneal;
    double speedTolerance = cConstants->speed_threshold*currentAnneal;
    if(posTolerance < cConstants->posDominationTolerance){
        posTolerance = cConstants->posDominationTolerance;
    }
    if(speedTolerance < cConstants->speedDominationTolerance){
        speedTolerance = cConstants->speedDominationTolerance;
    }

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
            //CHecks for a valid error status, so duplicates can be reset to valid without worry of overriding other error statuses
            if(duplicateCheck(tempAllAdults[i], tempAllAdults[j], cConstants, currentAnneal) && tempAllAdults[j].errorStatus == VALID) {
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
    
    //Is true if A is at least equally as good as B for all objectives
    bool AisEqual = false;
    //Is true if A is better than B for at least one objective
    bool AisBetter = false;
    //tolerances used to determine the range of values considered equal
    //these are both currently set to 1e-14 AU, I don't think these need to be modified 
    //this tolerance is about 0.0015m, and I don't think we can go lower?

    double posTolerance = 5;
    double speedTolerance = 5;
    //double posTolerance = 0.0;
    //double speedTolerance = 0.0;
    //true if A posDiff "equals" B posDiff
    bool APosEqualsB = false;
    //true if A speedDiff "equals" B speedDiff
    bool ASpeedEqualsB = false;

    //std::cout << "\n\nPERSON A TEST: posDiff = " << personA.posDiff << ", speedDiff = " << personA.speedDiff; 
    //std::cout << "\n\nPERSON B TEST: posDiff = " << personB.posDiff << ", speedDiff = " << personB.speedDiff; 

    //True is A posdiff is equal to B posDiff +- posTolerance
    if ((personA.posDiff < personB.posDiff + posTolerance) && (personA.posDiff > personB.posDiff - posTolerance)){
        APosEqualsB = true;
    }
    //True is A speeddiff is equal to B speedDiff +- speedTolerance
    if ((personA.speedDiff < personB.speedDiff + speedTolerance) && (personA.speedDiff > personB.speedDiff - speedTolerance)){
        ASpeedEqualsB = true;
    }
    //If the mission type is a rendezvous, A's speed diff needs to be lower for it to be equal or better
    if (missionType == Impact) {
        /*//Check if A's posDiff is approximately/better than B's posDiff and check the same for speed
        //If so, A is equal to B
        if ((personA.posDiff < personB.posDiff || APosEqualsB) && (personA.speedDiff > personB.speedDiff || ASpeedEqualsB)) {
            AisEqual = true; 
        }  
        //If A has a better posDiff or speedDiff than B, then A is better than B
        if (personA.posDiff < personB.posDiff || personA.speedDiff > personB.speedDiff) {
            AisBetter = true;
        }*/
        //For impact, only optimizing for one opbjective
        if ((personA.posDiff < personB.posDiff || APosEqualsB)) {
            AisEqual = true; 
        }  
        //If A has a better posDiff or speedDiff than B, then A is better than B
        if (personA.posDiff < personB.posDiff) {
            AisBetter = true;
        }
    }
    //If the mission type is a impact, A's speed diff needs to be higher than B's for it to be equal or better
    else {
        //If A.posDiff is approximately/better than B.posDiff, and A.speedDiff is approximately/better than B.speedDiff, then A is equal to B.
        if ((personA.posDiff < personB.posDiff || APosEqualsB) && (personA.speedDiff < personB.speedDiff || ASpeedEqualsB)) {
            AisEqual = true;
        }
        //If A has a better posDiff or speedDiff than B, then A is better than B
        if (personA.posDiff < personB.posDiff || personA.speedDiff < personB.speedDiff){
            AisBetter = true;
        }
        
    }

    //A Dominates B
    //A only dominates B if:
        //A.posDiff - posTolerance < B.posDiff <= A.posDiff, and A.speedDiff < B.speedDiff
        //A.posDiff < B.posDiff, and A.speedDiff - speedTolerance < B.speedDiff <= A.speedDiff
        //A.posDiff < B.posDiff and A.speedDiff < B.speedDiff (better in every way)
    if (AisEqual && AisBetter){
        //std::cout << "\n\nDomination check returns TRUE";
        return true;
    }
    //A and B are codominant if:
        //A.posDiff < B.posDiff & A.speedDiff > B.speedDiff
        //A.posDiff > B.posDiff & A.speedDiff < B.speedDiff
        //A.posDiff = B.posDiff & A.speedDiff = B.speedDiff (Unlikely since they're doubles, unless A & B are clones or genetically similar)
    //B dominates A if any of the conditions for "A dominates B", when reversed, are true.
    else {
        //std::cout << "\n\nDomination check returns FALSE";
        return false;
    }
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