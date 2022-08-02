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
  - It also is not practical to com