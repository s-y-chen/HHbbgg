#!/bin/sh

#Print out all bash commands
set -x

#Abort bash script on any error
set -e

#Print some basic debugging info
echo "whoami="`whoami`
echo "pwd="`pwd`
echo "hostname="`hostname`
echo "date="`date`
env

NJOB=$1
dataset=$2
year=$3

read -r workdir < path.txt

cd $workdir

export SCRAM_ARCH=slc7_amd64_gcc700
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`

eval `cmsenv`
env

ulimit -s unlimited

ref="data"

echo "to find raw ntuples for AOD ntuple starting from " $NJOB for $dataset "files in total"
if [ $dataset = $ref ]
then
    echo "analyzeHHbbgg processing data"
    ./analyzeHHbbgg Job${NJOB}_list.txt ${dataset}Job${NJOB}.root $dataset T T $year 
else
    echo "analyzeHHbbgg processing mc"
    ./analyzeHHbbgg Job${NJOB}_list.txt ${dataset}Job${NJOB}.root $dataset F T $year
fi 
