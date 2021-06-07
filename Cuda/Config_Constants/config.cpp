#include <iostream> // For cout
#include <fstream> // For file reading
#include <string>
#include <iomanip> // For setprecision in << operator
#include <time.h> // for time(0)
#include <math.h> // for sqrt() in constructor to derive v_escape from c3energy
#include "constants.h" // for AU

// Constructors uses geneticFileRead() to set the struct's properties from a default config file located in same folder as executable
cudaConstants::cudaConstants() {
    // Get values from the file
    FileRead("genetic.config");
    // Now that dry_mass and fuel_mass have been acquired, derive wet_mass
    this->wet_mass = this->dry_mass + this->fuel_mass;
    // Now that c3scale and c3energy have been assigned values, derive the final c3energy and v_escape
    this->c3energy *= this->c3scale;
    this->v_escape = sqrt(this->c3energy)/AU;
    // Assign cpu_numsteps to be equivalent to max_numsteps
    this->cpu_numsteps = this->max_numsteps;
}

// Operates same as default, however uses configFile as address for where the config file to be used is located
cudaConstants::cudaConstants(std::string configFile) {
    // Get values from the file
    FileRead(configFile);
    // Now that dry_mass and fuel_mass have been acquired, derive wet_mass
    this->wet_mass = this->dry_mass + this->fuel_mass;
    // Now that c3scale and c3energy have been assigned values, derive the final c3energy and v_escape
    this->c3energy *= this->c3scale;
    this->v_escape = sqrt(this->c3energy)/AU;
    // Assign cpu_numsteps to be equivalent to max_numsteps
    this->cpu_numsteps = this->max_numsteps;
}

// http://www.cplusplus.com/forum/beginner/11304/ for refesher on reading line by line
// Input: fileName - string address to the path to open the file being used
// Output - variable names found in the file that correspond to the cudaConstants' will be assigned the value followed by the name and '=', following a certain format assumption base (refer to config_readme.md for more precise information on format)
void cudaConstants::FileRead(std::string fileName) {
    // Use string line to hold a line read from the config file in variable configFile
    std::string line;
    std::ifstream configFile;
    configFile.open(fileName);

    if (configFile.is_open()) {
        // Go through line by line
        while ( std::getline(configFile, line) ) {
            // If line is not empty and the line is not staring as a comment, then the line is expected to be a variable constant being assigned 
            if (line != "" && ( line.find("//") != 0 )) {
                // Locate the "=" and use as pivot where prior is the variable name and after is variable value
                int equals_pivot = line.find("=");
                // Assumption made that there are no spaces from variable name to variable value, rather only occurring after a variable is assigned a value and afterwards may be an in-line comment
                int end_point = line.find_first_of(" ");

                // With the two positions acquired, capture the varable's name and the value it is being assigned
                std::string variableName = line.substr(0, equals_pivot   );
                std::string variableValue = line.substr( equals_pivot + 1, end_point - equals_pivot - 1);
                // Assign variableValue to the appropriate variable based on variableName, with conversion to the right data type
                // cudaConstant properties that are not expected in config are wet_mass, v_escape, and cpu_numsteps (those are to be derived in the constructor after this function is complete)
                if (variableName == "pos_threshold") {
                    this->pos_threshold = std::stod(variableValue);
                }
                else if (variableName == "coast_threshold") {
                    this->coast_threshold = std::stod(variableValue);
                }
                else if (variableName == "max_generations") {
                    this->max_generations = std::stoi(variableValue);
                }
                else if (variableName == "run_count") {
                    this->run_count = std::stoi(variableValue);
                }
                else if (variableName == "thruster_type") {
                    this->thruster_type = std::stoi(variableValue);
                }
                else if (variableName == "random_start") {
                    if (variableValue == "false") {
                        this->random_start = false;
                    }
                    else {
                        // If not set to false, then it is assumed the value is for true
                        this->random_start = true;
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
                else if (variableName == "record_mode") {
                    if (variableValue == "true") {
                        this->record_mode = true;
                    }
                    else {
                        // If not set to true, then it is assumed the value is false
                        this->record_mode = false;
                    }
                }
                else if (variableName == "initial_start_file_address") {
                    // Assumption that the address does not need to be converted/checked
                    this->initial_start_file_address = variableValue;
                }
                else if (variableName == "anneal_factor") {
                    this->anneal_factor = std::stod(variableValue);
                }
                else if (variableName == "write_freq") {
                    this->write_freq = std::stoi(variableValue);
                }
                else if (variableName == "disp_freq") {
                    this->disp_freq = std::stoi(variableValue);
                }
                else if (variableName == "best_count") {
                    this->best_count = std::stoi(variableValue);
                }
                else if (variableName == "change_check") {
                    this->change_check = std::stoi(variableValue);
                }
                else if (variableName == "anneal_initial") {
                    this->anneal_initial = std::stod(variableValue);
                }
                else if (variableName == "survivorRatio") {
                    this->survivorRatio = std::stod(variableValue);
                }
                else if (variableName == "mutation_rate") {
                    this->mutation_rate = std::stod(variableValue);
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
                else if (variableName == "coast_threshold") {
                    this->coast_threshold = std::stod(variableValue);
                }
                else if (variableName == "c3scale") {
                    this->c3scale = std::stod(variableValue);
                }
                else if (variableName == "c3energy") {
                    // Initially have c3energy just be the assigned value in the config, c3scale impacts c3energy (and by extension v_escape) within the constructor  
                    this->c3energy = std::stod(variableValue);
                }
                else if (variableName == "r_fin_ast") {
                    this->r_fin_ast = std::stod(variableValue);
                }
                else if (variableName == "theta_fin_ast") {
                    this->theta_fin_ast = std::stod(variableValue);
                }
                else if (variableName == "z_fin_ast") {
                    this->z_fin_ast = std::stod(variableValue);
                }
                else if (variableName == "vr_fin_ast") {
                    this->vr_fin_ast = std::stod(variableValue);
                }
                else if (variableName == "vtheta_fin_ast") {
                    this->vtheta_fin_ast = std::stod(variableValue);
                }
                else if (variableName == "vz_fin_ast") {
                    this->vz_fin_ast = std::stod(variableValue);
                }
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
                else if (variableName == "dry_mass") {
                   this->dry_mass = std::stod(variableValue);
                }
                else if (variableName == "fuel_mass") {
                    this->fuel_mass = std::stod(variableValue);
                }
                else if (variableName == "v_impact") {
                    this->v_impact = std::stod(variableValue);
                }
                else if (variableName == "rk_tol") {
                    this->rk_tol = std::stod(variableValue);
                }
                else if (variableName == "doublePrecThresh") {
                    this->doublePrecThresh = std::stod(variableValue);
                }
                else if (variableName == "GuessMaxPossibleSteps") {
                    this->GuessMaxPossibleSteps = std::stoi(variableValue);
                }
                else if (variableName == "max_numsteps") {
                    this->max_numsteps = std::stoi(variableValue);
                }
                else if (variableName == "min_numsteps") {
                    this->min_numsteps = std::stoi(variableValue);
                }
                else if (variableName == "num_individuals") {
                    this->num_individuals = std::stoi(variableValue);
                }
                else if (variableName == "survivor_count") {
                    this->survivor_count = std::stoi(variableValue);
                }
                else if (variableName == "thread_block_size") {
                    this->thread_block_size = std::stoi(variableValue);
                }
                else if (variableName == "timeRes") {
                    this->timeRes = std::stoi(variableValue);
                }
                else if (variableName == "time_seed") { // If the conifguration sets time_seed to NONE then time_seed is set to time(0) 
                    if (variableValue != "NONE") {
                        // If variableValue is not NONE, assumption is that it is a valid double value that can be converted and used
                        this->time_seed = std::stod(variableValue);
                    }
                    else {
                        this->time_seed = time(0);
                        std::cout << "time_seed value set to time(0)\n";
                    }
                }
                else {
                    // If none of the if cases were matches, then this is some unknown variable in the config file and output this to the terminal
                    std::cout << "Unknown variable '" << variableName <<"' in " << fileName <<"!\n";
                }
            }
        }
    }
    else {
        std::cout << "Unable to open " << fileName << " file!\n";
    }
}

// Output cudaConstant contents with formatting for better readibility when doing a run in main()
std::ostream& operator<<(std::ostream& os, const cudaConstants& object) {
    os << std::setprecision(12);
    os << "==========CONFIG=DATA===============================================================================\n";
    os << "Genetic Algorithm Related Values:\n";
    os << "\ttime_seed: "       << object.time_seed       << "\trandom_start: "   << object.random_start   << "\t\tnon_r_start_address: " << object.initial_start_file_address << "\n";
    os << "\tanneal_factor: "   << object.anneal_factor   << "\tanneal_initial: " << object.anneal_initial <<   "\tchange_check: "        << object.change_check << "\n";
    os << "\tnum_individuals: " << object.num_individuals << "\tthread_block_size: "     << object.thread_block_size << "\n";
    os << "\tsurvivor_count: "  << object.survivor_count  << "\tsurvivorRatio: "    << object.survivorRatio;
    os << "\tbest_count: "      << object.best_count      << "\t\tmax_generations: "<< object.max_generations<< "\trun_count: " << object.run_count << "\n\n";

    os << "Runge-Kutta Related Values:\n";
    os << "\trk_tol: " << object.rk_tol << "\tdoublePrecThresh: " << object.doublePrecThresh << "\ttimeRes: " << object.timeRes << "\tGuessMaxPossibleSteps: " << object.GuessMaxPossibleSteps << "\n";
    os << "\tmax_numsteps: " << object.max_numsteps << "\tmin_numsteps: "  << object.min_numsteps << "\tcpu_numsteps: " << object.cpu_numsteps << "\n\n";

    os << "Output Variables:\n";
    os << "\trecord_mode: " << object.record_mode << "\twrite_freq: " << object.write_freq << "\tdisp_freq: " << object.disp_freq << "\n\n";

    os << "Random Start Range Values:\n";
    os << "\tgamma: "  << object.gamma_random_start_range << "\ttau: " << object.tau_random_start_range << "\tcoast: " << object.coast_random_start_range << "\n";
    os << "\ttriptime min - max: " << object.triptime_min << " - "  << object.triptime_max << "\talpha: " << object.alpha_random_start_range << "\tbeta: " << object.beta_random_start_range << "\tzeta: " << object.zeta_random_start_range << "\n\n";
    
    os << "Mutation & Scales:\n";
    os << "\tmutation_rate: " << object.mutation_rate << "\n";
    os << "\tgamma_scale: "   << object.gamma_mutate_scale    << "\ttau_m_scale: "   << object.tau_mutate_scale   << "\tcoast_m_scale: " << object.coast_mutate_scale << "\n";
    os << "\talpha_m_scale: " << object.alpha_mutate_scale << "\tbeta_m_scale: "  << object.beta_mutate_scale  << "\tzeta_m_scale: " << object.zeta_mutate_scale << "\ttriptime_m_scale: "<< object.triptime_mutate_scale << "\n\n";

    os << "Spacecraft Info:\n";
    os << "\tthruster_type: " << object.thruster_type << "\tdry_mass: " << object.dry_mass << "\t\tfuel_mass: " << object.fuel_mass << "\t\twet_mass: " << object.wet_mass << "\n";
    // Display the c3energy assignment, showing c3scale to help clarify that it is not directly from config
    os << "\tc3energy (" << (object.c3scale * 100)    << "%): " << object.c3energy      << "\tv_escape: " << object.v_escape << "\t\tv_impact: " << object.v_impact << "\n";
    os << "\tpos_threshold: " << object.pos_threshold << "\tcoast_threshold: "<< object.coast_threshold<< "\n\n";

    os << "Asteriod Info:\n";
    os << "\t R: " << object.r_fin_ast  << "\t 0: " << object.theta_fin_ast << "\t Z: " << object.z_fin_ast << "\n";
    os << "\tvR: " << object.vr_fin_ast << "\tv0: " << object.vtheta_fin_ast << "\tvZ: " << object.vz_fin_ast << "\n\n";

    os << "Earth Info:\n";
    os << "\t R: " << object.r_fin_earth  << "\t 0: " << object.theta_fin_earth  << "\t Z: " << object.z_fin_earth  << "\n";
    os << "\tvR: " << object.vr_fin_earth << "\tv0: " << object.vtheta_fin_earth << "\tvZ: " << object.vz_fin_earth << "\n";
    os << "====================================================================================================\n";
    
    return os;
}
