add to sig_fit.C

#include "RooCrystalBall.h"

and inside sig_fit() function, at the beginning add 

gSystem->Load("RooCrystalBall_cxx.so");

//reference https://root.cern/manual/interacting_with_shared_libraries/
the to run, type on the terminal:

root -l

#note: if you have RooCrystalBall_cxx.so already, do not need to do the following step, go to gSystem->Load("RooCrystalBall_cxx.so") directly
.L RooCrystalBall.cxx++

gSystem->Load("RooCrystalBall_cxx.so")

.x sig_fit.C
