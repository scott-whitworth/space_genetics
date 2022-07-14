//Function that will adjust the annneal based on the previous anneal and if there was a change in the best individual
void changeAnneal (const std::vector<Adult>& oldAdults, const cudaConstants* cConstants, double & currentAnneal, const int & generation){
    
    //Grab the current best (based on rank-distance) progress of this generation
    double curProgress = oldAdults[0].progress;
    //this is used to make sure the anneal does not jump back up 
    double previousAnneal = currentAnneal;

    //This method changes the anneal based on the current progress
    //if the anneal reaches a new range, the anneal will step to a new value
    //This method currently has the fastest and most consistent convergence rate
    //Once it reaches 1e-3, the anneal will be function based on a power curve
    if(curProgress < 1e-8){ 
        currentAnneal = cConstants->anneal_initial; 
    } 
    else if(curProgress >= 1e-8 && curProgress < 1e-7){ 
        currentAnneal = cConstants->anneal_initial*.4; 
    } 
    else if(curProgress >= 1e-7 && curProgress < 1e-6){ 
        currentAnneal = cConstants->anneal_initial*.2; 
    } 
    else if(curProgress >= 1e-6 && curProgress < 1e-5){ 
        currentAnneal = cConstants->anneal_initial*.1; 
    } 
    else if(curProgress >= 1e-5 && curProgress < 1e-4){ 
        currentAnneal = cConstants->anneal_initial*.04; 
    } 
    else if(curProgress >= 1e-4 && curProgress < 1e-3){ 
        currentAnneal = cConstants->anneal_initial*.01; 
    }
    else{ 
        currentAnneal = (cConstants->anneal_initial * pow((1 - pow((curProgress), 0.07)) ,5)); 
        if(currentAnneal > previousAnneal){ 
            currentAnneal = previousAnneal; 
        } 
    }

    //Modify the current anneal by a sinusoidal multiplier
    //  This will allow some variability to occur within the anneal
    //  So, if the simulation is stuck, this will hopefully change things enough to move the simulation along
    //  Sin function modifies the current anneal ranging from +/- 20% of the current anneal, with the the percentage being based on the generation
    currentAnneal += currentAnneal * .7 * sin(generation * (M_PI/25));
     
    //Check to make sure that the current anneal does not fall below the designated minimum amount
    if (currentAnneal < cConstants->anneal_final)
    {
        currentAnneal = cConstants->anneal_final;
    }
}