root -l

#note: if you have RooCrystalBall_cxx.so already, do not need to do the following step, go to gSystem->Load("RooCrystalBall_cxx.so") directly
.L RooCrystalBall.cxx++

gSystem->Load("RooCrystalBall_cxx.so")

.x dimass_fit_to_ws.C

python make_cards.py
