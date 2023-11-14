#ifndef OBJECTIVE_H
#define OBJECTIVE_H

#include <string>

//enum used to store the potential parameter goals
//  *Important Note For Future*:
//      Negative values indicate that the program will try to minimize a value and
//      positive values indicate a favor for large values
enum parameterGoals {
    //Minimization goals
    POS_DIFF = 1, 
    SPEED_DIFF = 2,
    FUEL_SPENT = 3,
    TRIP_TIME = 4,
    ORBIT_POS_DIFF = 5,
    ORBIT_SPEED_DIFF = 6,
    MARS_DIST = 7, 
    HORZ_ANGLE_DIFF = 8,
    VERT_ANGLE_DIFF = 9,
    ORBIT_ASST = 10,

    //Error value
    UNKNOWN = 0,
};

//See mission_config_readme.md for more details

//This structure holds the necessary information for a program objective
//This allows for the program to handle dynamic program objectives

//NOTE: When adding a new objective, make sure to complete these steps:
//  1) Make sure the parameter is accessible from within the child class
//  2) Add the goal to the parameterGoals enum, making sure that if its a minimization objective, the set number is negative and vice versa
//  3) Add the new objective goal to the importObjective function within config.cpp
//  5) Add to child's getParameter function so that it returns the correct parameter
//  6) Add a new sort function to the adult.cpp to allow the program to sort by the objective's parameter in the desired order
//  7) Add to the parameterSort() function in sort.cpp so that the program uses the new sort when the function is called with the new goal

struct objective {
    //The name the program will use for the objective
    std::string name;

    //The parameter goal of the objective
    parameterGoals goal; 

    //The convergence threshold; this will be used to determine whether a parameter has been solved by the best individual
    // double convergenceThreshold;
    //The requested target
    double target; 

    //The maximum allowed difference, smaller difference will be considered solved for the objective
    double allowedDifference;
    //The difference that will be optimized for, the algorithm will not try and optimize past this value for the mission
    //  This determines when the program will/will not prioritize better values
    //      If two individuals compared in a domination check have parameter values better than the goal difference,
    //          the parameter will not be considered for domination
    double goalDifference; 

    //equateTolerance determines the differentiation needed for two parameters to be considered not equal
    //  this is primarily a floating point comparison issue
    double equateTolerance;

    //Default constructor, will set the goal as an error
    objective(): goal(UNKNOWN){}

    //Constructor that accepts the information for an objective, should be the only constructor used
    //    Should only be called when the mission.config file is parsed/loaded
    objective(std::string _name, parameterGoals _goal, double _target, double _allowedDiff, double _goalDiff, double _equateTolerance):
        name(_name), goal(_goal), target(_target), allowedDifference(_allowedDiff), goalDifference(_goalDiff), equateTolerance(_equateTolerance) {}
};

#endif