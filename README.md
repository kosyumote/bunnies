## Space-discrete stochastic predator-prey model

Here, all the source code for Konstantin Borisov and Tatyana Yatsenko's final project for Dr. Rubin's Mathematical Biology class can be found.

*Note that this project only works on Linux systems, due to some "system" calls, and the use of a Linux binary in the folder utils. However, the source code may be freely used to create an executable for other OSes.*

##### bunnies.cpp

This is the basic version of the model, which creates animated pictures of the predator and prey populations over time.

##### bunnies_initial_value.cpp

This version of the model creates a picture that indicates what state the model will end in with various initial conditions.

##### bunnies_diverge.cpp

This version of the model simulates the same initial conditions several times, recording the results in data files. *Note that this version of the model should be completely cross-platform, as it doesn't rely on trickery to achieve the desired result.*