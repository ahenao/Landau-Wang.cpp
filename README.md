# Landau-Wang.cpp
/*
This Program perform a single spin flip Glauber Dynamics simulation of the Ising
Magnetic system, by the technique of Monte Carlo

Andres Henao Aristizabal
Physics Engineer- Universidad Nacional de Colombia
Studying Master in Computational Physics- Universitat Politecnica Catalunya, Unversitat Barcelona
ahenaoa@unal.edu.co 
*/


The Program Landau-Wang.cpp is the main, do

g++ Landau-Wang.cpp mersenne.o userintf.o 
./a.out -L 12 

You can change 12 by the linear lattice size.
Then in screen will appear the option to the interval that you want. (Note: Inside the program, you can change the size of the interval, in screen you only choose the beggining of the interval. The program is made with a delta of energy of 0.01, so if you begin in the first interval -2.03 and the lenght you put it to meas_int=217 you will cover up to 0.13).

At the end, the output g1 will have the natural logarithm of the density of states in the second column and the energy in the first one.

