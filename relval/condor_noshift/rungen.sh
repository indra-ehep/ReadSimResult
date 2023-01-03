#!/bin/bash
#To be run on remote machine
#Take input arguments as an array
myArray=( "$@" )
#Array: Size=$#, an element=$1, all element = $@

printf "Start Running Histogramming at ";/bin/date
printf "Worker node hostname ";/bin/hostname
export HOME=/afs/cern.ch/user/i/idas
CMSVER=CMSSW_13_0_X_2022-12-31-1100
SCRAM_ARCH=slc7_amd64_gcc11

if [ -z ${_CONDOR_SCRATCH_DIR} ] ; then 
    echo "Running Interactively" ; 
else
    echo "Running In Batch"
    echo ${_CONDOR_SCRATCH_DIR}
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    scramv1 project CMSSW $CMSVER
    cd $CMSVER/src
    eval `scramv1 runtime -sh`
    #git cms-addpkg Geometry/HGCalCommonData
    git cms-merge-topic 40404
    scram b -j 8
    #git cms-addpkg Configuration/Generator
    #cd ../..
    cp ../../generator.tar.gz .
fi

pwd
ls -la
echo "Unzipping tar zipped package"
tar --strip-components=0 -zxvf generator.tar.gz

if [ -z ${_CONDOR_SCRATCH_DIR} ] ; then 
    echo "Running Interactively" ; 
else
    cd relval
    echo "Current working directory : " ${PWD}
    pwd
    ls -la
    echo "Before copy"
    ls -la ../Geometry/HGCalCommonData/src/HGCalDDDConstants.cc
    cp HGCalDDDConstants.cc ../Geometry/HGCalCommonData/src/HGCalDDDConstants.cc
    echo "After copy"
    ls -la ../Geometry/HGCalCommonData/src/HGCalDDDConstants.cc
    scram b 
fi
#Run for Base, Signal region

echo "All arguements: "$@
echo "Number of arguements: "$#
geom=$1
cshift=$2
index=$3

echo "Before copy"
ls -la ../Geometry/HGCalCommonData/data/hgcalHEmix/v17/hgcalHEmix.xml
#cp hgcalHEmix_cs_${cshift}.xml ../Geometry/HGCalCommonData/data/hgcalHEmix/v17/hgcalHEmix.xml
echo "After copy"
ls -la ../Geometry/HGCalCommonData/data/hgcalHEmix/v17/hgcalHEmix.xml
echo "Has changed"

cmsRun testHGCalSIMSingleMuonPt100_cfg.py geometry=$geom
mv step1.root step1_${index}.root
mv step1.log step1_${index}.log

ls -ltr

pwd

printf "Simulation completed at ";/bin/date
#---------------------------------------------
#Copy the ouput root files
#---------------------------------------------
condorOutDir1=/eos/cms/store/group/dpg_hgcal/comm_hgcal/geomval/cassette_noshift_v17_20230103/$geom/$cshift
#condorOutDir1=/eos/user/i/idas/SimOut/geomval/muontomo_newtrig/$geom
condorOutDir=/cms/store/user/idas/SimOut/geomval/cassette_noshift_v17_20230103/$geom/$cshift

if [ -z ${_CONDOR_SCRATCH_DIR} ] ; then
    echo "Running Interactively" ;
else
    xrdcp -f step1_${index}.root root://eosuser.cern.ch/${condorOutDir1}
    xrdcp -f step*_${index}.log root://eosuser.cern.ch/${condorOutDir1}
    xrdcp -f step*_${index}.* root://se01.indiacms.res.in:1094/${condorOutDir}
    echo "Cleanup"
    cd ../../
    rm -rf $CMSVER
    rm *.root
fi
printf "Done ";/bin/date

