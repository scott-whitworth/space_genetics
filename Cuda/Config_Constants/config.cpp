#include <iostream> // For cout
#include <fstream> // For file reading
#include <string>
#include <iomanip> // For setprecision in << operator
#include <time.h> // for time(0)
#include <math.h> // for sqrt() in constructor to derive v_escape from c3energy
#include <cctype> // for tolower(), needed for Tesla machine only

// Default constructor which assumes that genetic.config and mission.config are the files being pulled from
// Asteroid file is determined within configFile
cudaConstants::cudaConstants() {
    //Set initial algorithm state to unspecified
    //  If the user specifies an algorithm, it will be set from the genetic.config file read
    algorithm = UNSPECIFIED;

    // Get values from the genetic config file
    FileRead("../Config_Constants/genetic.config");
    // Get values from the mission config file
    FileRead("../Config_Constants/mission.config"); 

    //get the destination
    FileRead("../Config_Constants/" + this->destination);

    //Check if there was no algorithm specified
    if (algorithm == UNSPECIFIED) {
        //Set the algorithm based on the number of objectives
        if (missionObjectives.size() > 2) {
            //Rank rarity is best for missions with many objectives
            std::cout << "\nAutomatically setting algorithm to rank-rarity.\n";
            algorithm = RANK_RARITY;
        }
        else {
            //Rank rarity is best for missions with few objectives
            std::cout << "\nAutomatically setting algorithm to rank-distance.\n";
            algorithm = RANK_DISTANCE;
        }
    }

    //If time mutation scale is not set, set it to the difference between triptime max and min
    if (this->triptime_mutate_scale < this->doublePrecThresh) {
        //std::cout << "\nTEST: setting triptime mutate scale to difference\n";
        this->triptime_mutate_scale = (this->triptime_max - this->triptime_min);
    }

    //If time mutation scale is not set, set it to the difference between triptime max and min
    if (this->carryover_individuals > this->num_individuals) {
        this->carryover_individuals = 0;
    }

    // Now that dry_mass and fuel_mass have been acquired, derive wet_mass
    this->wet_mass = this->dry_mass + this->fuel_mass;

    // Now that c3scale and c3energy have been assigned values, derive the final c3energy and v_escape
    this->c3energy *= this->c3scale;
    this->v_escape = sqrt(this->c3energy)/AU;
}



// http://www.cplusplus.com/forum/beginner/11304/ for refesher on reading line by line
// Input: fileName - string address to the path to open the file being used
// Output - variable names found in the file that correspond to the cudaConstants' will be assigned the value followed by the 
// name and '=', following a certain format assumption base (refer to config_readme.md for more precise information on format)
void cudaConstants::FileRead(std::string fileName) {
    // Use string line to hold a line read from the config file in variable configFile
    std::string line;
    std::ifstream configFile;
    configFile.open(fileName);

    if (configFile.is_open()) {
        // Go through line by line
        while ( std::getline(configFile, line) ) {
            // If line is not empty and the line is not starting as a comment, then the line is expected to be a variable constant being assigned 
            if (line != "" && ( line.find("//") != 0 )) {

                //First check to see if mission objectives are being imported, since that will be a different process than other config variables
                if (line == "Mission_Objectives:")
                {
                    //Get the next line to get the first objective
                    std::getline(configFile, line); 

                    //Import objectives until the section is done
                    while (line != "" && line.find("//") != 0 ) {
                        //Call the import objective function to import the info from the objective line
                        importObjective(line); 

                        //Get the next objective
                        std::getline(configFile, line); 
                    }
                }
                //Import a normal config file variable
                else {
                    // Locate the "=" and use as pivot where prior is the variable name and after is variable value
                    int equals_pivot = line.find("=");
                    // Assumption made that there are no spaces from variable name to variable value, rather only occurring after a variable is assigned a value and afterwards may be an in-line comment
                    int end_point = line.find_first_of(" ");

                    // With the two positions acquired, capture the varable's name and the value it is being assigned
                    std::string variableName = line.substr(0, equals_pivot);
                    std::string variableValue = line.substr( equals_pivot + 1, end_point - equals_pivot - 1);
                    // Assign variableValue to the appropriate variable based on variableName, with conversion to the right data type
                    // cudaConstant properties that are not expected in config are wet_mass, v_escape (those are to be derived in the constructor after this function is complete)

                    //Convert the variable value to lower case for easier processing
                    for (int i = 0; i < variableValue.size(); i++) {
                        variableValue[i] = std::tolower(variableValue[i]); 
                    }


//////////////////////////////////////////////////////////////////// -- INITIALIZING & RANDOM -- /////////////////////////////////////////////////////////////////////////
                    if (variableName == "time_seed") { // If the conifguration sets time_seed to NONE then time_seed is set to time(0) 
                        if (variableValue != "none") {
                            // If variableValue is not NONE, assumption is that it is a valid double value that can be converted and used
                            this->time_seed = std::stod(variableValue);
                        }
                        else {
                            this->time_seed = time(0);
                            std::cout << "time_seed value set to time(0)\n";
                        }
                    }
                    else if (variableName == "max_generations") {
                        this->max_generations = std::stoi(variableValue);
                    }
                    else if (variableName == "run_count") {
                        this->run_count = std::stoi(variableValue);
                    }
                    else if (variableName == "carryover_individuals") {
                        this->carryover_individuals = std::stoi(variableValue);
                    }


//////////////////////////////////////////////////////////////////////// -- RUNGE KUTTA -- /////////////////////////////////////////////////////////////////////////////
                    else if (variableName == "rk_tol") {
                        this->rk_tol = std::stod(variableValue);
                    }
                    else if (variableName == "doublePrecThresh") {
                        this->doublePrecThresh = std::stod(variableValue);
                    }
                    else if (variableName == "max_numsteps") {
                        this->max_numsteps = std::stoi(variableValue);
                    }
                    else if (variableName == "min_numsteps") {
                        this->min_numsteps = std::stoi(variableValue);
                    }


//////////////////////////////////////////////////////////////////// -- POOL & THREAD BLOCK -- /////////////////////////////////////////////////////////////////////////
                    else if (variableName == "num_individuals") {
                        this->num_individuals = std::stoi(variableValue);
                    }
                    else if (variableName == "survivor_count") {
                        this->survivor_count = std::stoi(variableValue);
                    }
                    else if (variableName == "thread_block_size") {
                        this->thread_block_size = std::stoi(variableValue);
                    }


////////////////////////////////////////////////////////////////////////// -- OUTPUT -- /////////////////////////////////////////////////////////////////////////////////
                    else if (variableName == "record_mode") {
                        if (variableValue == "true") {
                            this->record_mode = true;
                        }
                        else {
                            // If not set to true, then it is assumed the value is false
                            this->record_mode = false;
                        }
                    }

                    else if (variableName == "write_freq") {
                        this->write_freq = std::stoi(variableValue);
                    }
                    else if (variableName == "all_write_freq") {
                        this->all_write_freq = std::stoi(variableValue);
                    }
                    else if (variableName == "disp_freq") {
                        this->disp_freq = std::stoi(variableValue);
                    }

                
///////////////////////////////////////////////////////////////////////// -- SPACECRAFT -- //////////////////////////////////////////////////////////////////////////////
                    else if (variableName == "thruster_type") {
                        this->thruster_type = std::stoi(variableValue);
                    }
                    else if (variableName == "dry_mass") {
                        this->dry_mass = std::stod(variableValue);
                    }
                    else if (variableName == "fuel_mass") {
                        this->fuel_mass = std::stod(variableValue);
                    }


/////////////////////////////////////////////////////////////////////////// -- MISSION -- ////////////////////////////////////////////////////////////////////////////////
                    else if (variableName == "destination") {
                        this->destination = variableValue;
                    }
                    else if (variableName == "c3scale") {
                        this->c3scale = std::stod(variableValue);
                    }
                    else if (variableName == "c3energy") {
                        // Initially have c3energy just be the assigned value in the config, c3scale impacts c3energy (and by extension v_escape) within the constructor  
                        this->c3energy = std::stod(variableValue);
                    }


////////////////////////////////////////////////////////////////////// -- GENETIC ALGORITHM -- ///////////////////////////////////////////////////////////////////////////
                    else if (variableName == "algorithm_type") {
                        //Assign the algorithm type based on the user input
                        if (variableValue == "rank-rarity") {
                            this->algorithm = RANK_RARITY;
                            std::cout << "\nRank-Rarity algorithm specified\n";
                        }
                        else if (variableValue == "rank-distance") {
                            this->algorithm = RANK_DISTANCE;
                            std::cout << "\nRank-Distance algorithm specified\n";
                        }
                        else {
                            std::cout << "\nAlgorithm selection unidentified\n";
                        }
                    }
                    else if (variableName == "alpha_random_start_range") {
                        this->alpha_random_start_range = std::stod(variableValue);
                    }
                    else if (variableName == "beta_random_start_range") {
                        this->beta_random_start_range = std::stod(variableValue);
                    }
                    else if (variableName == "zeta_random_start_range") {
                        this->zeta_random_start_range = std::stod(variableValue);
                    }
                    else if (variableName == "triptime_max") {
                        this->triptime_max = std::stod(variableValue) * SECONDS_IN_YEAR;
                    }
                    else if (variableName == "triptime_min") {
                        this->triptime_min = std::stod(variableValue) * SECONDS_IN_YEAR;
                    }
                    else if (variableName == "gamma_random_start_range") {
                        this->gamma_random_start_range = std::stod(variableValue);
                    }
                    else if (variableName == "tau_random_start_range") {
                        this->tau_random_start_range = std::stod(variableValue);
                    }
                    else if (variableName == "coast_random_start_range") {
                        this->coast_random_start_range = std::stod(variableValue);
                    }

                    else if (variableName == "mutation_amplitude") {
                        this->mutation_amplitude = std::stod(variableValue);
                    }

                    else if (variableName == "default_mutation_chance") {
                        this->default_mutation_chance = std::stod(variableValue);
                    }


                    else if (variableName == "gamma_mutate_scale") {
                        this->gamma_mutate_scale = std::stod(variableValue);
                    }
                    else if (variableName == "tau_mutate_scale") {
                        this->tau_mutate_scale = std::stod(variableValue);
                    }
                    else if (variableName == "coast_mutate_scale") {
                        this->coast_mutate_scale = std::stod(variableValue);
                    }
                    else if (variableName == "triptime_mutate_scale") {
                        this->triptime_mutate_scale = std::stod(variableValue) * SECONDS_IN_YEAR;
                    }
                    else if (variableName == "zeta_mutate_scale") {
                        this->zeta_mutate_scale = std::stod(variableValue);
                    }
                    else if (variableName == "beta_mutate_scale") {
                        this->beta_mutate_scale = std::stod(variableValue);
                    }
                    else if (variableName == "alpha_mutate_scale") {
                        this->alpha_mutate_scale = std::stod(variableValue);
                    }
                    else if (variableName == "anneal_initial") {
                        this->anneal_initial = std::stod(variableValue);
                    }
                    else if (variableName == "anneal_final") {
                        this->anneal_final = std::stod(variableValue);
                    }

/////////////////////////////////////////////////////////////////////////// -- OTHER -- //////////////////////////////////////////////////////////////////////////////////
                    else if (variableName == "timeRes") {
                        this->timeRes = std::stoi(variableValue);
                    }
                    else if (variableName == "sun_r_min") {
                        this->sun_r_min = std::stod(variableValue);
                    }
                    else if (variableName == "maxSimNum") {
                        this->maxSimNum = std::stoi(variableValue);
                    }
                    else if (variableName == "best_count") {
                        this->best_count = std::stoi(variableValue);
                    }
                    else if (variableName == "orbitalRadius"){
                        this->orbitalRadius = std::stod(variableValue);
                    }
                    else if (variableName == "orbitalSpeed"){
                        this->orbitalSpeed = std::stod(variableValue);
                    }
                    else if (variableName == "MSOI_scale"){
                        this->MSOI_scale = std::stod(variableValue);
                    }
                    else if (variableName == "coast_threshold") {
                        this->coast_threshold = std::stod(variableValue);
                    }
                    else if (variableName == "divisions") {
                        this->divisions = std::stoi(variableValue);
                    }


/////////////////////////////////////////////////////////////////////////// -- TARGET -- //////////////////////////////////////////////////////////////////////////////////
                    else if (variableName == "r_fin_target") {
                        this->r_fin_target = std::stod(variableValue);
                    }
                    else if (variableName == "theta_fin_target") {
                        this->theta_fin_target = std::stod(variableValue);
                    }
                    else if (variableName == "z_fin_target") {
                        this->z_fin_target = std::stod(variableValue);
                    }
                    else if (variableName == "vr_fin_target") {
                        this->vr_fin_target = std::stod(variableValue);
                    }
                    else if (variableName == "vtheta_fin_target") {
                        this->vtheta_fin_target = std::stod(variableValue);
                    }
                    else if (variableName == "vz_fin_target") {
                        this->vz_fin_target = std::stod(variableValue);
                    }
                    else if (variableName == "orbitalPeriod") {
                        this->orbitalPeriod = stod(variableValue);
                    }
                    else if (variableName == "gravAssistDist") {
                        this->gravAssistDist = stod(variableValue);
                    }
                    else if (variableName == "gravAssistTime") {
                        this->gravAssistTime = stod(variableValue);
                    }
                    else if (variableName == "gravAssistTimeFrac") {
                        this->gravAssistTimeFrac = stod(variableValue);
                    }


/////////////////////////////////////////////////////////////////////////// -- EARTH -- //////////////////////////////////////////////////////////////////////////////////
                    else if (variableName == "r_fin_earth") {
                        this->r_fin_earth = std::stod(variableValue);
                    }
                    else if (variableName == "theta_fin_earth") {
                        this->theta_fin_earth = std::stod(variableValue);
                    }
                    else if (variableName == "z_fin_earth") {
                        this->z_fin_earth = std::stod(variableValue);
                    }
                    else if (variableName == "vr_fin_earth") {
                        this->vr_fin_earth = std::stod(variableValue);
                    }
                    else if (variableName == "vtheta_fin_earth") {
                        this->vtheta_fin_earth = std::stod(variableValue);
                    }
                    else if (variableName == "vz_fin_earth") {
                        this->vz_fin_earth = std::stod(variableValue);
                    }

/////////////////////////////////////////////////////////////////////////// -- MARS -- //////////////////////////////////////////////////////////////////////////////////
                    else if (variableName == "r_fin_mars") {
                        this->r_fin_mars = std::stod(variableValue);
                    }
                    else if (variableName == "theta_fin_mars") {
                        this->theta_fin_mars = std::stod(variableValue);
                    }
                    else if (variableName == "z_fin_mars") {
                        this->z_fin_mars = std::stod(variableValue);
                    }
                    else if (variableName == "vr_fin_mars") {
                        this->vr_fin_mars = std::stod(variableValue);
                    }
                    else if (variableName == "vtheta_fin_mars") {
                        this->vtheta_fin_mars = std::stod(variableValue);
                    }
                    else if (variableName == "vz_fin_mars") {
                        this->vz_fin_mars = std::stod(variableValue);
                    }

/////////////////////////////////////////////////////////////////////////// -- ERROR -- //////////////////////////////////////////////////////////////////////////////////
                    else {
                        // If none of the if cases were matches, then this is some unknown variable in the config file and output this to the terminal
                        std::cout << "Unknown variable '" << variableName <<"' in " << fileName <<"!\n";
                    }
                }
            }
        }
    }
    else {
        std::cout << "Unable to open " << fileName << " file!\n";
    }
}

// Transfers information from a config file line to the mission objectives vector
void cudaConstants::importObjective(std::string line) {
    
    //Create two pivot points to get information between the comma-based boundaries
    int beginningPivot, endPivot; 

    //temp storage variables that will be used to create the new objective object
    std::string name;
    parameterGoals goal; 
    double target, diff, equateTolerance; 

    //Temp string will assist with eliminating spaces from the line and identifying the goal
    std::string tempStr; 


    //Remove any spaces from the line using the temp string
    for (int i = 0; i < line.size(); i++) {
        if (line[i] != ' ') {
            tempStr.push_back(line[i]);
        }
    }
    //temp string now equals line without the spaces, make equate line to temp string
    line = tempStr;

    //Find the first end pivot to get the name
    endPivot = line.find(",");

    //Pull the name from the substring of the line
    name = line.substr(0, endPivot);

    //Convert the string to lower case for easier decernment of parameter goals
    //Name has already been taken, so no issue if it is changed in the line
    for (int i = 0; i < line.size(); i++) {
        line[i] = std::tolower(line[i]); 
    }

    //Get the next pivot points
    beginningPivot = endPivot+1;
    endPivot = line.find(",", beginningPivot);

    //temp string used for the parameter goal
    tempStr = line.substr(beginningPivot, endPivot - beginningPivot); 

    //Determine the parameter goal based on the next segment of the imported line
    if (tempStr == "min_pos_diff") {
        //Will optimize for minimum position difference
        goal = MIN_POS_DIFF;
    }
    else if (tempStr == "min_speed_diff") {
        //optimize for the minimum speed difference
        goal = MIN_SPEED_DIFF; 
    }
    else if (tempStr == "min_horz_velocity_diff") {
        //optimize for a final horizontal velocity angle difference 
        goal = MIN_HORZ_VEL_DIFF; 
    }
    else if (tempStr == "min_vert_velocity_diff") {
        //optimize for a final vertical velocity angle difference 
        goal = MIN_VERT_VEL_DIFF; 
    }
    else if (tempStr == "min_fuel_spent") {
        //optimize for minimal fuel usage
        goal = MIN_FUEL_SPENT;
    }
    else if (tempStr == "min_trip_time") {
        //Optimize for minimal trip time
        goal = MIN_TRIP_TIME;
    }
    else if (tempStr == "min_orbit_pos_diff") {
        goal = MIN_ORBIT_POS_DIFF;
    }
    else if (tempStr == "min_orbit_speed_diff") {
        goal = MIN_ORBIT_SPEED_DIFF;
    }
    else if (tempStr == "min_mars_dist") {
        goal = MIN_MARS_DIST;
    }
    else if (tempStr == "max_speed_diff") {
        //Optimize for highest speed
        goal = MAX_SPEED_DIFF; 
    }
    else if (tempStr == "max_orbit_asst"){
        //Optimize for the highest change in angular momentum during an assist
        goal = MAX_ORBIT_ASST;
    }
    else {
        //No parameter goal identified
        goal = UNKNOWN; 
    }

    //Find the next pivot point for the convergence tolerance
    beginningPivot = endPivot+1;
    endPivot = line.find(",", beginningPivot); 

    //Pull the convergence tolerance from the substring
    target = std::stod(line.substr(beginningPivot, endPivot - beginningPivot + 1));

    //Find the next pivot point for the domination tolerance
    beginningPivot = endPivot+1;
    endPivot = line.find(",", beginningPivot);
    
    //Get the domination threshold from the last substring
    diff = std::stod(line.substr(beginningPivot, endPivot - beginningPivot + 1));

    //Find the last pivot points
    beginningPivot = endPivot+1; 
    endPivot = line.size(); 

    //Pull the equateTolerance
    equateTolerance = std::stod(line.substr(beginningPivot, endPivot - beginningPivot + 1));

    //Checks for the user to make sure they entered the goal, the convergence threshold, and the domination threshold correctly
    //See if the goal is set correctly
    if (goal == 0) {
        std::cout << "\n-----BAD OBJECTIVE GOAL PULLED; BAD OBJECTIVE: " << name << "-----\n";
    }
    // //Check to see if the goal is to minimize, but the domination threshold is set higher than the convergence threshold
    // else if (goal < 0 && dominationThreshold > convergenceThreshold) {
    //     std::cout << "\n-----DOMINATION THRESHOLD SET TOO HIGH; BAD OBJECTIVE: " << name << "-----\n";
    // }
    // //Check to see if the goal is to maximize, but the domination threshold is set lower than the convergence threshold
    // else if (goal > 0 && dominationThreshold < convergenceThreshold) {
    //     std::cout << "\n-----DOMINATION THRESHOLD SET TOO LOW; BAD OBJECTIVE: " << name << "-----\n";
    // }

    //See if the objective is a maximization
    //If so, set the thresholds as a 1/thresholds so the rest of the code can treat it as a minimization
    // if (goal > 0) {
    //     convergenceThreshold = 1/convergenceThreshold;
    //     dominationThreshold = 1/dominationThreshold;
    // }    

    //Add the objective to the mission objectives vector using the gathered information
    missionObjectives.push_back(objective(name, goal, target, diff, equateTolerance)); 
}

// Output cudaConstant contents with formatting for better readibility when doing a run in main()
std::ostream& operator<<(std::ostream& os, const cudaConstants& object) {
    os << std::setprecision(12);
    os << "\n==========CONFIG=DATA===============================================================================\n";
    os << "Genetic Algorithm Related Values:\n";
    os << "\ttime_seed: "       << object.time_seed       << "\tnum_divisions: "    << object.divisions         << "\talgorithm: "          << object.algorithm         << "\n";
    os << "\tnum_individuals: " << object.num_individuals << "\tsurvivor_count: "   << object.survivor_count    << "\tthread_block_size: "  << object.thread_block_size << "\n";
    os << "\tbest_count: "      << object.best_count      << "\tmax_generations: "  << object.max_generations   << "\trun_count: "          << object.run_count         << "\n\n";

    os << "Runge-Kutta Related Values:\n";
    os << "\trk_tol: "       << object.rk_tol       << "\tdoublePrecThresh: " << object.doublePrecThresh << "\ttimeRes: "   << object.timeRes    << "\n";
    os << "\tmax_numsteps: " << object.max_numsteps << "\tmin_numsteps: "     << object.min_numsteps     << "\tsun_r_min: " << object.sun_r_min  << "\n\n";

    os << "Output Variables:\n";
    os << "\trecord_mode: " << object.record_mode << "\twrite_freq: " << object.write_freq << "\tdisp_freq: " << object.disp_freq << "\n\n";

    os << "Random Start Range Values:\n";
    os << "\tgamma: "        << object.gamma_random_start_range << "\ttau: "  << object.tau_random_start_range  << "\tcoast: " << object.coast_random_start_range << "\n";
    os << "\talpha: "        << object.alpha_random_start_range << "\tbeta: " << object.beta_random_start_range << "\tzeta: "  << object.zeta_random_start_range  << "\n";
    os << "\ttriptime min: " << object.triptime_min             << "\tmax: "  << object.triptime_max            << "\n\n";
    
    os << "Mutation & Scales:\n";
    os << "\tmutation_amplitude: " << object.mutation_amplitude     << "\tanneal_initial: " << object.anneal_initial    << "\n";
    os << "\tgamma_scale: "        << object.gamma_mutate_scale     << "\ttau_m_scale: "    << object.tau_mutate_scale  << "\tcoast_m_scale: " << object.coast_mutate_scale << "\n";
    os << "\talpha_m_scale: "      << object.alpha_mutate_scale     << "\tbeta_m_scale: "   << object.beta_mutate_scale << "\tzeta_m_scale: "  << object.zeta_mutate_scale  << "\n";
    os << "\ttriptime_m_scale: "   << object.triptime_mutate_scale  << "\n\n";

    os << "Spacecraft Info:\n";
    os << "\tdry_mass: " << object.dry_mass << "\t\tfuel_mass: " << object.fuel_mass << "\t\twet_mass: " << object.wet_mass << "\n";
    // Display the c3energy assignment, showing c3scale to help clarify that it is not directly from config
    os << "\tc3energy ("        << (object.c3scale * 100) << "%): "         << object.c3energy << "\tv_escape: "    << object.v_escape  << "\n";
    os << "\tthruster_type: "   << object.thruster_type   << "\tcoast_threshold: " << object.coast_threshold << "\n\n";

    os << "Target Info:\n";
    os << "\t R: " << object.r_fin_target  << "\t 0: " << object.theta_fin_target  << "\t Z: " << object.z_fin_target  << "\n";
    os << "\tvR: " << object.vr_fin_target << "\tv0: " << object.vtheta_fin_target << "\tvZ: " << object.vz_fin_target << "\n\n";

    os << "Earth Info:\n";
    os << "\t R: " << object.r_fin_earth  << "\t 0: " << object.theta_fin_earth  << "\t Z: " << object.z_fin_earth  << "\n";
    os << "\tvR: " << object.vr_fin_earth << "\tv0: " << object.vtheta_fin_earth << "\tvZ: " << object.vz_fin_earth << "\n\n";

    os << "Mars Info:\n";
    os << "\t R: "         << object.r_fin_mars  << "\t 0: " << object.theta_fin_mars  << "\t Z: " << object.z_fin_mars  << "\n";
    os << "\tvR: "         << object.vr_fin_mars << "\tv0: " << object.vtheta_fin_mars << "\tvZ: " << object.vz_fin_mars << "\n";
    os << "\tMSOI_scale: " << object.MSOI_scale  << "\n\n";

    os << "Mission Goals: ";
    for (int i = 0; i < object.missionObjectives.size(); i++) {
        os << "\n\tObjective: " << object.missionObjectives[i].name << "\tIdentified Goal: " << object.missionObjectives[i].goal << "\tTarget: " << object.missionObjectives[i].target 
           << "\tAllowed Difference: " << object.missionObjectives[i].allowedDifference << "\tEquate Tolerance: " << object.missionObjectives[i].equateTolerance; 
    }
    os << "\n";
    os << "====================================================================================================\n";
    
    return os;
}
