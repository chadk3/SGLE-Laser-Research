# SGLE-Laser-Research
Some of my graduate research involving models and search algorithms for high energy mode locked laser systems based off of the Sinusoidal Ginzburg Landau Equation. Includes simulation software to model the output of a given set of parameters, a test for mode locking, a genetic algorithm to seek mode locked extrema, and a toroidal search algorithm to map out the 5+ dimensional dynamics with an adaptive time step to assure convergence. There are also some basic grid search files that enabled me to find local extrema to use as initial conditions for the genetic algorithm, and some programs that plot the transmission curves.

Sinusoidal_Ginzburg_Landau_Equation.m is the simulation program. It relies on sgle_rhs.m

ga_sgle.m is the genetic algorithm

Toroidal_Sampling_SGLE.m is the toroidal search algorithm with an adaptive time step.
