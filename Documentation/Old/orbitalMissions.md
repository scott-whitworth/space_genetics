<h1>Current progress of orbital missions </h1>

<h2>Psyche mission </h2>
Currently this mission does not work at all and there are many factors as to why that may be. One major factor that the orbital mission and the Odyssey mission aim to examine is the problem of resolution when it comes to time and stepSize. The code is likely taking too large of steps for a mission like this that uses the gravity of Mars. These values of 1600 for max_numsteps, anbd 400 for min_numsteps need to be larger, and 3600 for timeRes should be smaller when including Mars' gravity. timeRes should be around 450-900 which will make the planet calculations take longer but it is necessary. The max_numsteps should be as large as possible (3000 for the computers in the lab) and min_numsteps should be at least 800. These will make the overall calculations more accurate but not the calculations in the MSOI. To make these more accurate, MSOI_steps and MSOI_error will likely need to be implemented and tested (or callRk will need to be called twice). Another reason this mission may not work is potentially an algorithmic problem, however this is unlikely (it may be worth looking into). The last and most problematic reason has to do with something currently out of our control. The last time we ran the psyche mission, the spacecraft would orbit the sun for the entire triptime and would not use mars. The very few times it did use mars it would end up going far too fast. If this is still the case, it means that the thrusters we are using are not capable of completing a mission like this (this is more likely for an orbital mission as this one has more time to slowdown). 

<h2> Orbital missions (Mars) </h2>

<h3> marsOrbit.config </h3>

<h4>Background</h4>
The marsOrbit config uses the expected flyby date for the psyche mission with Mars. This mission was made for a gravitational assist and not orbit, so it will likely never work. The few times this mission was tested, there was a lot of error in the velocities and calculation of the conservation of mechanical energy. Like the psyche mission, the spacecraft is unable to slow down enough or even attempt to orbit mars. The expectation of the results currently should be that the mission will achieve the posDiff, but not the speedDiff. The speedDiff will be too fast because of Mars' gravity and it is unknown on how to fix this currently.

<h4>Tests conducted</h4>
One of the things that have been tested regarding this mission is adding more steps within MSOI to make sure the gravitational calculations are accurate. The stepSize could easily step past important calculations for this mission at the end of it, so its important to have as many steps as possible at the end. Through the tests done, this seems to help the error go down at least but does not improve speedDiff. Another test that has been done is changing the c3 and its scale to see if had too much or not enough launch energy. The values for c3 are not based off anything specific for this mission, so they should be changed to some extent. Lastly, the timeRes was made very small (225) and this also led to no major change other than making the code slow. 

<h3> marsOdyssey.config </h3>

<h4>Background</h4>
This config is based off the Odyssey mission conducted by NASA in 2001. This mission was sent to Mars to orbit around it (and study it). The initial orbit it would maintain was much larger than the final orbit it was meant to maintain (the largest orbit had a period of ~20 hours while the final one was ~2 hours). This was done by using the spacecraft's powerful thrusters (much more powerful than the ones in our code) to readjust its orbital period. The Odyssey's closest approach was ~330km above the surface of Mars. The final orbit it maintained was a circular orbit, but the other orbits were all very elliptical. With all of this in mind, this mission may prove difficult for our code to simulate. 

<h4>Tests conducted</h4>

- The major change to better this mission was changing the way vtheta_Rate was calculated in motion_equations.cpp. When calculating the gravity of Mars, it is important to use the fmod of the two thetas being used. Initially, it was being calculated as "fmod(y.theta, 2pi) - fmod(mars.theta, 2pi)", but it should be "fmod(fmod(y.theta, 2pi) - fmod(mars.theta, 2pi))". This changed the outcome for this mission greatly when it came to the posDiff and speedDiff, but did not affect the marsOrbit mission. The change was important, however it led to the posDiff getting stuck at around 0.001, which was not happening before. Without this change, the posDiff will converge, but the speedDiff will be around 1e-6, and the tangential velocity spikes and so does the error. 
- To get the posDiff unstuck, more data points are needed, which is where the MSOI_steps and MSOI_error come in. On my personal computer with a very powerful GPU, we were able to do runs that had 20,000+ data points. The extra data points came from making MSOI_steps be 100 and MSOI_error be 5.0. By doing this the posDiff would slowly get down to as low as 0.0004 (only 400 generations were ran). The major takeaway from this run was not the posDiff but the error and tangential velocity. By making the MSOI much larger, the error was in th 1e-3 range and the tangential velocity was finally smooth. These graphs can be seen on teams in a document called MSOI_graphs in results->mars_orbit_results.
- Another important change made for these tests was removing the MAX_DISTANCE given to any individuals that met a domination tolerance. Including this seemingly hindered progress. However this was only an initial observation so it may be worth messing with. 
- Within output.cpp in printFinalGen(), make sure to include the bestProgress sort done before finalRecord() so the individual that is used for analysis is not the best speed individual but the best position individual.
- c3 scale was generally kept at 1.0, but its worth messing with. 
- Some tests were conducted where the mutate scale and random start were changed for zeta. We tried making it zero for both or only positive (for the Odyssey mission it would be better if it was only negative), and this did not help much. 

<h4>Future changes</h4>

- To implement more points on a machine with a less powerful GPU, callRk will likely have to be called twice. The initial call, and then again when it reaches some factor of MSOI. By doing this the watchdog problem can be avoided and more points can be calculated.
- After implementing more points, it is important to actually see the orbit. Currently, the trajectory ends when the closest individual to the convergence threshold are met for posDiff and speedDiff. It is important to actually see that the spacecraft is orbiting Mars as we hope. 
- Different thrusters (maybe). This may only be relevant for the psyche mission (it may not be), but its something Dr. K briefly looked into for the psyche mission at least.
- Double check all of the motion equations that deal with Mars.

<h3> Setting up an Orbit mission </h3>

- In the mission.config set destination=marsOdyssey.config or destination=marsOrbit.config, and use the values for the orbitPosDiff and orbitSpeedDiff
- The triptimes for an orbit mission should not exceed 0.5-1.0, it takes the average mission 7 months to complete (ours should take longer because they might need a very slow approach)
- The only things to adjust in the two orbit configs are MSOI_error and MSOI_steps (c3 can be adjusted too). The MSOI_error should probably not exceed 2.0 (depending on the max_numsteps and device). As for MSOI_steps, if MSOI_error is 2.0, then 20 is probably appropriate (depending on the max_numsteps and device). If MSOI_error is 1.0, 200 is a good value for MSOI_steps (the maxs steps is roughly 3100 for the lab computers, so keep an eye on bestCount)
- In the genetic.config, make the min_numsteps 800, timeRes 900 (or lower), max_numsteps should probably stay as 1600 unless MSOI_error or MSOI_steps are small. Mutation scale for triptime should match the difference of max_triptime - min_triptime (for all missions). 
- Make sure the bestProgress individual is being selected at the end of each run
- It may be a good idea to comment out the section in sort.cpp that starts with: for (int j = 1; j < validAdults-1; j++) {. This section gives all the individuals under the convergence threshold distance + 1, this does not help for this mission (probably). 
- If the code crashes after many generations, it is likely due to the callRk taking too many steps and MSOI_error or MSOI_steps needs to be lowered. 