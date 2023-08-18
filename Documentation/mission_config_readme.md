<h1>Mission.config objective Explanations</h1>

<h2>Overview</h2>

- Mission.config allows you to set objectives dynamically
- The genetic algorithim will attempt to optimize for each set objective

<br>

# Setting Objectives

- In mission.config, objectives are set in the lines after "Mission_Objectives:"
- Each line consists of one objective
  - An objective consists of a name, the goal, the convergence threshold, the domination threshold, and the equate tolerance
  - The structure of an objective is: *name*, *goal*, *convergence threshold*, *domination threshold*, *equate tolerance*
- The line after the last objective must be empty or have a comment

<br>

<h2>Objective Name</h2>

- The name of the objective, which is primarily used for data reporting
- There are no rules when setting the objective name, however it is valuable to make it match the objective goal for clarity

<br>

<h2>Objective Goal</h2>

- The objective goal tells the genetic algorithim what output to optimize for
- The current list of goals are:
  - <b>Min_Pos_Diff</b> (minimizes final position difference between the indivdual and the target)
  - <b>Min_Speed_Diff</b> (minimized final speed difference between the indivdual and the target)
  - <b>Max_Speed_Diff</b> (maximizes final speed difference between the indivdual and the target)
  - <b>Min_Fuel_Spent</b> (minimizes the fuel spent during the simulation)
  - <b>Min_Trip_Time</b> (minimizes the trip time of the simulation)
  - <b>Min_Orbit_Pos_Diff</b> (minimizes the final position difference between the individual and the orbit radius of the target body)
  - <b>Min_Orbit_Speed_Diff</b> (minimizes the final speed difference between the individual and the orbit speed of the target body)
- Note: Goals which optimize for input variables (e.g. min_trip_time) currently don't work nearly as well and may be better optimized by setting challenging min/max values
- Note: Orbit goals don't currently work very well, see [orbitalMissions.md](orbitalMissions.md) for detail

<br>

<h2>Convergence Tolerance</h2>

- The convergence tolerance determines the value which the objective is considered solved
- If an individual has variables better than all objectives' convergence tresholds, it is considered converged
- If the objective goal is to minimize an output, the objective is considered solved when the individual's variable is less than the objective's convergence threshold
- If the objective goal is to maximize an output, the objective is considered solved when the individual's variable is greater than the objective's convergence threshold

<br>

<h2>Domination Threshold</h2>

- The domination tolerance is the value which determines if two individuals' variables are compared for the objective
- When two individuals are compared to determine domination, if both individuals' variables are better than the domination threshold, the objective is not factored into the comparison
- Similarly to convergence threshold, if the goal is a minimization, the two individuals' variables have to be smaller than the domination threshold for the objective not to be considered for domination, and vice versa for maximization

<br>

<h2>Equate Tolerance</h2>

- This value determines how similar two individuals' variables have to be for the individuals to be considered in the objective
- Used because the objectives have different units
  - It also is not practical to compare parameters at unreasonable precisions (e.g. measuring fuel spent at 1e-14 kg)

<br>
</br>

# Building New Objectives

- The current structure of how objectives are handled makes the process of building new objectives relatively simple
- Candidates for new objectives ideally should be a value that can be minimized or maximized and do not directly have to do with a single input parameter
- Following are the steps needed to build a new objective:
    1. Initialize the output variable that will be used to measure the progress the objective within the [child class](../Genetic_Algorithm/child.h)
       - How this variable is calculated may be from within the runge-kutta simulation (see fuelSpent) or calculated within its own function (see posDiff)
    2. Add the in-code name of the goal to the parameterGoals enum in [objective.h](objective.h)
       - Important Note: if this objective is about minimizing an output, make sure the value of the enum variable is <0 (and vice-versa for maximizing objectives). The sign of objective enum variable is a quick method used by the code to determine the desired direction of the objective. 
    3. Using the new enum goal, add the new objective to the section in the importObjective() function in [config.cpp](config.cpp) where the inputs from mission.config are translated to objective objects.
    4. Add the new objective to the getParameter() function in child.cpp
       - This function is used to dynamically return the correct child class variable depending on the objective
       - Thus, you need to add an if statement which returns the output variable you created in step 1 if the input objective is the one you are building
    5. Build a new sort function for the objective's output variable in the adult class
       - These sort functions are used when calculating an individual's distance
       - The sort should return true when the first adult submitted has a better value for the objectives output variable than the second adult
         - For example, if the objective is a minimization, the sort function needs to favor smaller values
         - See previously created sort functions and their associated objectives for clarity
    6. Finally, add the sort that you just built to the parameterSort() function in sort.cpp
       - This function calls the correct sort based on the submitted objective
       - So, append to the switch statement within the function to sort based on the new objective using the sort function from the previous step 