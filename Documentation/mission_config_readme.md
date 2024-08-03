<h1>Mission.config objective Explanations</h1>

<h2>Overview</h2>

- Mission.config allows you to set objectives without recompiling the code
- The genetic algorithim will attempt to optimize for each set objective

<br>

# Setting Objectives

- In mission.config, objectives are set in the lines after "Mission_Objectives:"
- Each line consists of one objective
  - An objective consists of a goal, the target value, the allowed difference (convergence threshold), the goal difference (domination threshold), and the equate tolerance
  - The structure of an objective is: *goal*, *target value*, *allowed difference*, *goal difference*, *equate tolerance*
- The line after the last objective must be empty or have a comment

<br>


<h2>Objective Goal</h2>

- The objective goal tells the genetic algorithim what results to optimize for
- The current list of goals are:
  - <b>Pos_Diff</b> (Target final position difference between the target and the craft [Au] )
  - <b>Speed_Diff</b> (Target final speed difference between the target and the craft [Au/s] )
  - <b>Orbit_Pos_Diff</b> (Target final difference between the orbital distance from the target and the craft [Au] )
  - <b>Orbit_Speed_Diff</b> (Target final difference between the desired orbital velocity around the target and the craft's velocity [Au/s] )
  - <b>Horz_Angle_Diff</b> (Target final horizontal plane velocity angle difference between the target and the craft [degrees] )
  - <b>Vert_Angle_Diff</b> (Target final vertical plane velocity angle difference between the target and the craft [degrees] )
  - <b>Fuel_Spent</b> (Target fuel used by the spacecraft during the simulation [Kg] )
  - <b>Trip_Time</b> (Target trip time of the simulations [s] )

<br>

<h2>Target Value</h2>

- The target value tells the genetic algorithm what goal value to optimize towards

<br>

<h2>Allowed Difference</h2>

- The allowed difference is the (+/-) range around the target value for which the objective is considered solved
- If an individual has objective differences better than all objectives' convergence tresholds, it is considered converged


<br>

<h2>Goal Difference</h2>

- The goal difference (AKA domination tolerance) is the (+/-) range which determines when to stop optimizing for that specific objective.
- When two individuals are compared to determine domination, if both individuals' objective differences are better (lower) than the domination threshold, the objective is not factored into the comparison


<br>

<h2>Equate Tolerance</h2>

- This value determines how precise the values have to be to be considered equal
- Used because the objectives have different units
  - It also is not practical to compare parameters at unreasonable precisions (e.g. measuring fuel spent at 1e-14 kg)

<br>
</br>
