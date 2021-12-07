#!/bin/bash

#To be run on remote machine
#Take input arguments as an array
myArray=( "$@" )
#Array: Size=$#, an element=$1, all element = $@

printf "Start Running Histogramming at ";/bin/date
printf "Worker node hostname ";/bin/hostname
CMSVER=CMSSW_12_1_1

if [ -z ${_CONDOR_SCRATCH_DIR} ] ; then 
    echo "Running Interactively" ; 
else
    echo "Running In Batch"
    echo ${_CONDOR_SCRATCH_DIR}
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    SCRAM_ARCH=slc7_amd64_gcc900
    scramv1 project CMSSW $CMSVER
    cd $CMSVER/src
    eval `scramv1 runtime -sh`
    mkdir -p Configuration
    cp -r $CMSSW_RELEASE_BASE/src/Configuration/Generator  Configuration/
    cp ../../SingleMuPt100_hgcal_cfi.py Configuration/Generator/python/SingleMuPt100_hgcal_cfi.py
fi

if [ -z ${_CONDOR_SCRATCH_DIR} ] ; then 
    echo "Running Interactively" ; 
else
    pwd
    ls -la
    scram b -j 4
fi
#Run for Base, Signal region

echo "All arguements: "$@
echo "Number of arguements: "$#
geom=$1
index=$2

cmsDriver.py SingleMuPt100_hgcal_cfi -s GEN,SIM -n 10000 --conditions auto:phase2_realistic_T21 --beamspot HGCALCloseBy --datatier GEN-SIM --eventcontent FEVTDEBUG --geometry $geom --era Phase2C11I13M9 --relval 9000,100 --fileout file:step1_${index}.root  --customise_commands process.RandomNumberGeneratorService.generator.initialSeed="cms.untracked.uint32($RANDOM)" --no_exec
grep -n "process.RandomNumberGeneratorService.generator.initialSeed" SingleMuPt100_hgcal_cfi_GEN_SIM.py
cmsRun SingleMuPt100_hgcal_cfi_GEN_SIM.py

ls -ltr

pwd

printf "Simulation completed at ";/bin/date
#---------------------------------------------
#Copy the ouput root files
#---------------------------------------------
condorOutDir=/eos/user/i/idas/MuDeltaPt100GeV
if [ -z ${_CONDOR_SCRATCH_DIR} ] ; then
    echo "Running Interactively" ;
else
    xrdcp -f step1_${index}.root root://eosuser.cern.ch/${condorOutDir}
    echo "Cleanup"
    cd ../../
    rm -rf $CMSVER
    rm *.root
fi
printf "Done ";/bin/date
