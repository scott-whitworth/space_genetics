<h1>Current progress of orbital missions </h1>

<h2>Psyche mission </h2>
Currently this mission does not work at all and there are many factors as to why that may be. One major factor that the orbital mission and the Odyssey mission aim 
to examine is the problem of resolution when it comes to time and stepSize. These are both likely far too large for a mission like this that uses the gravity of mars
at their current values. These values of 1600 for max_numsteps, 400 for min_numsteps, and 3600 for timeRes need to be much smaller when including Mars' gravity.
timeRes should be around 450 which will make the planet calculations take longer but it is neccessary. The max_numsteps should be as large as possible (3000 for the 
computers in the lab) and min_numsteps should be at least 800. These will make the overall calculations more accurate but not the calculations in the MSOI. To make 
these more accurate, MSOI_steps and MSOI_error will likely need to be implemented and tested. Another reason this mission may not work is potentially an algorithmic problem 
however this is unlikely (it may be worth looking into). The last and most problematic reason has to do with something currently out of our control. The last time we ran the psyche
mission, the spacecraft would orbit the sun for the entire triptime and would not use mars. The very few times it did use mars it would end up going far too fast. If this is still the case, it means that the thrusters we are using are not capable of completing a mission like this (this is more likely for an orbital mission as this one has more time to slow
down). 

<h2> Orbital missions (Mars) </h2>
