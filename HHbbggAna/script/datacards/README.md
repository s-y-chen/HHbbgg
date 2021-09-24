root -l

#note: if you have RooCrystalBall_cxx.so already, do not need to do the following step, go to gSystem->Load("RooCrystalBall_cxx.so") directly
.L RooCrystalBall.cxx++

gSystem->Load("RooCrystalBall_cxx.so")

.x dimass_fit_to_ws.C

python make_cards.py

####

script to run limit:

1. first set up Higgs combine:

https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/

following the instructions in the link above:

```
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_2_13
cd CMSSW_10_2_13/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
```

2. use the command to produce limit:

```
combine -M AsymptoticLimits -m 125 -n versionname HHbbgg_datacard.txt --saveWorkspace --saveToys --run blind
```
