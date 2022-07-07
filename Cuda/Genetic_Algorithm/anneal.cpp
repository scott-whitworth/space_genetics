//Function that will adjust the annneal based on the previous anneal and if there was a change in the best individual
void changeAnneal (const std::vector<Adult>& oldAdults, const cudaConstants* cConstants, double & currentAnneal, int & oldestBirthday, double & dRate, const int & generation, double & previousProgress){
    
    //Grab the current best (based on rank-distance) progress of this generation
    double curProgress = oldAdults[0].progress;
    //this is used to make sure the anneal does not jump back up 
    double previousAnneal = currentAnneal;
    

    //Caluclate the new current anneal
    //It will be adjusted based on a power scale

    //How we calculate the anneal depends on the mission type
    //  Our current formula is anneal_initial * (1 - curProgress^0.1) * dRate
    //  dRate is used to help move anneal along when there isn't a ton of progress being made
    //  The section in the parentheses should equal 1 when the goals has been met
    //  However, to make sure the anneal doesn't go too low, we make sure that it doen't go below the min anneal
    //  Thus, the number of mission goals are different between the missions, so we need to check what type of mission it is to calculate the anneal

    if (cConstants -> missionType == Impact) 
    {
        currentAnneal = cConstants->anneal_initial * (1 - (curProgress));
    }
    else 
    {   
        if(previousProgress > curProgress){
            //Check to see if the new progress is worse than previous generations
            //  This means that a new best rank-distance adult has been found with a lower progress, increasing the anneal to a figure that is too large
            //reset the progress so it does not jump anneals
            curProgress = previousProgress;  
        }else{//the progress is the same or better
            //update previousProgress to be the new progress
            previousProgress = curProgress;
        }

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
            currentAnneal = cConstants->anneal_initial*.4; 
        } 
        else if(curProgress >= 1e-6 && curProgress < 1e-5){ 
            currentAnneal = cConstants->anneal_initial*.09; 
        } 
        else if(curProgress >= 1e-5 && curProgress < 1e-4){ 
            currentAnneal = cConstants->anneal_initial*.04; 
        } 
        else if(curProgress >= 1e-4 && curProgress < 1e-3){ 
            currentAnneal = cConstants->anneal_initial*.008; 
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
        currentAnneal += currentAnneal * .2 * sin(generation * (M_PI/30));
        

        //OUTDATED - keeping here for reference
        //Here is the best anneal method that uses a sine function
        /*
        if(curProgress < 1e-8){ 
            currentAnneal = cConstants->anneal_initial; 
        } 
        else if(curProgress > 1e-8 && curProgress < 1e-6){ 
            currentAnneal = cConstants->anneal_initial*.3 + .1*cConstants->anneal_initial*sin(generation*(M_PI/25)); 
        } 
        else if(curProgress > 1e-6 && curProgress < 1e-5){ 
            currentAnneal = cConstants->anneal_initial*.2 + .1*cConstants->anneal_initial*sin(generation*(M_PI/25)); 
        } 
        else if(curProgress > 1e-5 && curProgress < 1e-4){ 
            currentAnneal = cConstants->anneal_initial*.1 + .08*cConstants->anneal_initial*sin(generation*(M_PI/25)); 
        } 
        else if(curProgress > 1e-4 && curProgress < 1e-3){ 
            currentAnneal = cConstants->anneal_initial*.05 + .025*cConstants->anneal_initial*sin(generation*(M_PI/25)); 
        } 
        else{ 
            currentAnneal = (cConstants->anneal_initial * pow((1 - pow((curProgress), 0.07)),5))*dRate; 
            if(currentAnneal > previousAnneal){ 
                currentAnneal = previousAnneal; 
            } 
        }
        */

    } 
    //Check to make sure that the current anneal does not fall below the designated minimum amount
    if (currentAnneal < cConstants->anneal_final)
    {
        currentAnneal = cConstants->anneal_final;
    }
}