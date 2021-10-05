#!/bin/bash
#To be run on remote machine
#Take input arguments as an array
myArray=( "$@" )
#Array: Size=$#, an element=$1, all element = $@

printf "Start Running Histogramming at ";/bin/date
printf "Worker node hostname ";/bin/hostname

if [ -z ${_CONDOR_SCRATCH_DIR} ] ; then 
    echo "Running Interactively" ; 
else
    echo "Running In Batch"
    echo ${_CONDOR_SCRATCH_DIR}
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    export SCRAM_ARCH=slc7_amd64_gcc900
    scramv1 project CMSSW CMSSW_12_0_2
    cd CMSSW_12_0_2/src
    eval `scramv1 runtime -sh`
    git cms-addpkg Configuration/Generator
    cd ../..
    
fi

tar --strip-components=1 -zxvf generator.tar.gz
cp GenConfig/SingleMuPt100_hgcal_cfi.py CMSSW_12_0_2/src/Configuration/Generator/python/SingleMuPt100_hgcal_cfi.py

if [ -z ${_CONDOR_SCRATCH_DIR} ] ; then 
    echo "Running Interactively" ; 
else
    cd CMSSW_12_0_2/src
    scram b -j 4
fi
#Run for Base, Signal region
#./complib.sh
echo "All arguements: "$@
echo "Number of arguements: "$#
geom=$1
index=$2

cmsDriver.py SingleMuPt100_hgcal_cfi -s GEN,SIM -n 10 --conditions auto:phase2_realistic_T21 --beamspot HGCALCloseBy --datatier GEN-SIM --eventcontent FEVTDEBUG --geometry $geom --era Phase2C11I13M9 --relval 9000,100 --fileout file:step1_${index}.root  --customise_commands process.RandomNumberGeneratorService.generator.initialSeed="cms.untracked.uint32($RANDOM)" --no_exec #--nThreads 4
grep -n "process.RandomNumberGeneratorService.generator.initialSeed" SingleMuPt100_hgcal_cfi_GEN_SIM.py
cmsRun SingleMuPt100_hgcal_cfi_GEN_SIM.py


printf "Simulation completed at ";/bin/date
#---------------------------------------------
#Copy the ouput root files
#---------------------------------------------
condorOutDir1=/eos/user/i/idas/SimOut/DeltaPt/$geom
condorOutDir=/cms/store/user/idas/SimOut/DeltaPt/$geom
if [ -z ${_CONDOR_SCRATCH_DIR} ] ; then
    echo "Running Interactively" ;
else
    #xrdcp -f ${sample}_tree_*.root root://se01.indiacms.res.in:1094/${condorOutDir}/${year} 
    xrdcp -f step1_${index}.root root://eosuser.cern.ch/${condorOutDir1}/${geom}
    echo "Cleanup"
    rm -rf CMSSW_10_2_14
    rm *.root
fi
printf "Done ";/bin/date
