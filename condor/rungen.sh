#!/bin/bash
#To be run on remote machine
#Take input arguments as an array
myArray=( "$@" )
#Array: Size=$#, an element=$1, all element = $@

printf "Start Running Histogramming at ";/bin/date
printf "Worker node hostname ";/bin/hostname
CMSVER=CMSSW_12_5_0_pre5
export HOME="/afs/cern.ch/user/i/idas"
echo "Home set as :"
echo ${HOME}

if [ -z ${_CONDOR_SCRATCH_DIR} ] ; then 
    echo "Running Interactively" ; 
else
    echo "Running In Batch"
    echo ${_CONDOR_SCRATCH_DIR}
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    SCRAM_ARCH=slc7_amd64_gcc10
    scramv1 project CMSSW $CMSVER
    cd $CMSVER/src
    eval `scramv1 runtime -sh`
    git cms-addpkg Configuration
    #cd ../..
    cp ../../generator.tar.gz .
fi

pwd
ls -la
tar --strip-components=0 -zxvf generator.tar.gz
cp ReadSimResult/GenConfig/SingleMuPt100_hgcal_cfi.py Configuration/Generator/python/SingleMuPt100_hgcal_cfi.py

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
nevent=10000

# cmsDriver.py SingleMuPt100_hgcal_cfi -s GEN,SIM -n 100000 --conditions auto:phase2_realistic_T21 --beamspot HGCALCloseBy --datatier GEN-SIM --eventcontent FEVTDEBUG --geometry $geom --era Phase2C11I13M9 --relval 9000,100 --fileout file:step1_${index}.root  --customise_commands process.RandomNumberGeneratorService.generator.initialSeed="cms.untracked.uint32($RANDOM)" --no_exec  --nThreads 4
# grep -n "process.RandomNumberGeneratorService.generator.initialSeed" SingleMuPt100_hgcal_cfi_GEN_SIM.py
# cmsRun SingleMuPt100_hgcal_cfi_GEN_SIM.py

cmsDriver.py SingleMuPt100_hgcal_cfi  -s GEN,SIM -n $nevent --conditions auto:phase2_realistic_T21 --beamspot HLLHC14TeV --datatier GEN-SIM --eventcontent FEVTDEBUG  --geometry $geom --era Phase2C11I13M9 --relval 9000,100 --fileout file:step1_${index}.root --customise_commands process.RandomNumberGeneratorService.generator.initialSeed="cms.untracked.uint32($RANDOM)" --no_exec --nThreads 8 
grep -n "process.RandomNumberGeneratorService.generator.initialSeed" SingleMuPt100_hgcal_cfi_GEN_SIM.py
cmsRun SingleMuPt100_hgcal_cfi_GEN_SIM.py > step1_${index}.log  2>&1


cmsDriver.py step2  -s DIGI:pdigi_valid,L1TrackTrigger,L1,DIGI2RAW,HLT:@fake2 --conditions auto:phase2_realistic_T21 --datatier GEN-SIM-DIGI-RAW -n $nevent --eventcontent FEVTDEBUGHLT --geometry $geom --era Phase2C11I13M9 --filein  file:step1_${index}.root  --fileout file:step2_${index}.root  --nThreads 8 > step2_${index}.log  2>&1
 
cmsDriver.py step3  -s RAW2DIGI,RECO,RECOSIM,PAT,VALIDATION:@phase2Validation+@miniAODValidation,DQM:@phase2+@miniAODDQM --conditions auto:phase2_realistic_T21 --datatier GEN-SIM-RECO,MINIAODSIM,DQMIO -n $nevent --eventcontent FEVTDEBUGHLT,MINIAODSIM,DQM --geometry $geom --era Phase2C11I13M9  --filein  file:step2_${index}.root  --fileout file:step3_${index}.root  --nThreads 8 > step3_${index}.log  2>&1


ls -ltr
pwd

printf "Simulation completed at ";/bin/date
#---------------------------------------------
#Copy the ouput root files
#---------------------------------------------
#condorOutDir1=/eos/cms/store/group/dpg_hgcal/comm_hgcal/geomval/etaphi_debug/$geom
condorOutDir1=/eos/user/i/idas/SimOut/geomval/etaphi_debug/$geom
condorOutDir=/cms/store/user/idas/SimOut/geomval/etaphi_debug/$geom
if [ -z ${_CONDOR_SCRATCH_DIR} ] ; then
    echo "Running Interactively" ;
else
    xrdcp -f step*_${index}.root root://eosuser.cern.ch/${condorOutDir1}
    xrdcp -f step*_${index}.log root://eosuser.cern.ch/${condorOutDir1}
    xrdcp -f step*_${index}.* root://se01.indiacms.res.in:1094/${condorOutDir}
    echo "Cleanup"
    cd ../../
    rm -rf $CMSVER
    rm *.root
fi
printf "Done ";/bin/date
