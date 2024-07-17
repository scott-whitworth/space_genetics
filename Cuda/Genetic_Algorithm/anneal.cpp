
double changeIndividualAnneal (double curProgress, const cudaConstants* cConstants, const int & generation){
    

    //this is used to make sure the anneal does not jump back up 
    double currentAnneal;
    //This method changes the anneal based on the current progress
    //if the anneal reaches a new range, the anneal will step to a new value
    //This method currently has the fastest and most consistent convergence rate
    //Once it reaches 1e-3, the anneal will be function based on a power curve
    
    //This function is here because it has the needed shape through the range 1e-3 - 1.0 
    currentAnneal = (cConstants->anneal_initial * pow((1 - pow((curProgress), 0.08)) ,5)); 
    //Do not want the anneal to go backwards at this point, make sure it only decreases


    //Modify the current anneal by a sinusoidal multiplier
    //  This will allow some variability to occur within the anneal
    //  So, if the simulation is stuck, this will hopefully change things enough to move the simulation along
    //  Sin function modifies the current anneal ranging from +/- 70% of the current anneal, with the the percentage being based on the generation
    //  This function has no mathematical significance other than it yeilds the best convergence rate
    currentAnneal += currentAnneal * .7 * sin(generation * (M_PI/25));
    // currentAnneal += currentAnneal * .05 * sin(generation * (M_PI/25));
     
    //Check to make sure that the current anneal does not fall below the designated minimum amount
    if (currentAnneal < cConstants->anneal_final)
    {
        currentAnneal = cConstants->anneal_final;
    }
    return currentAnneal;
}