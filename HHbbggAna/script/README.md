add to sig_fit.C

#include "RooCrystalBall.h"

and inside sig_fit() function, at the beginning add 

gSystem->Load("RooCrystalBall_cxx.so");

//reference https://root.cern/manual/interacting_with_shared_libraries/
the to run, type on the terminal:

root -l

.L RooCrystalBall.cxx++

gSystem->Load("RooCrystalBall_cxx.so")

.x sig_fit.C
