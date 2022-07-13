#ifndef OBJECTIVE_H
#define OBJECTIVE_H

#include <string>

//enum used to store the potential parameter goals
//  *Important Note For Future*
//   Negative values indicate that the program will try to minimize a value and positive values indicate a favor for large values
enum parameterGoals {
    //Minimization goals
    MIN_POS_DIFF = -1, 
    MIN_SPEED_DIFF = -2,
    MIN_FUEL_SPENT = -3,
    MIN_TRIP_TIME = -4,

    //Error value
    UNKNOWN = 0, 

    //Maximization goals
    MAX_SPEED_DIFF = 2, 
};

//This structure holds the necessary information for a program objective
//This will allow for the program to handle dynamic program objectives

//TODO: When this is complete, write out a quick guide on what to build where when implementing a new objective
//      Specifically, things will need to be added to the adult function that returns the right value and new sorts will need to be added
struct objective {
    //The name the program will use for the objective
    std::string name;

    //The parameter goal of the objective
    parameterGoals goal; 

    //The convergence threshold; this will be used to determine whether a parameter has been solved by the best individual
    double convergenceThreshold; 

    //The domination threshold determines when the program will/will not prioritize better values
    //      If two individuals compared in a domination check have parameter values better than the domination check, the parameter will not be considered for domination
    double dominationThreshold; 

    //equateTolerance determines the differentiation needed for two parameters to be considered not equal
    double equateTolerance;

    //Default constructor, will set the goal as an error
    objective(): goal(UNKNOWN){}

    //Constructor that accepts the information for an objective, should be the only constructor used
    objective(std::string _name, parameterGoals _goal, double _convergenceThreshold, double _dominationThreshold, double _equateTolerance):
        name(_name), goal(_goal), convergenceThreshold(_convergenceThreshold), dominationThreshold(_dominationThreshold), equateTolerance(_equateTolerance) {}
};

#endif