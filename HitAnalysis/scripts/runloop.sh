#!/bin/bash

#inputdir=/home/idas/t3store3/root_files/HGCAL_Geometry/SimOut/DeltaPt/Extended2026D86
#inputdir=/eos/user/i/idas/SimOut/DeltaPt/Extended2026D86
inputdir=/eos/user/p/psuryade/SimOut/DeltaPt/Extended2026D86_25
pydir=$PWD/ReadSimResult/SimTrackAna/python
for i in `seq 0 9`
do
  echo Processing loop $i with file step1_${i}.root
  # if [ -f $inputdir/step1.root ] ; then
  #     rm $inputdir/step1.root
  # fi
  #ln -s $inputdir/step1_${i}.root $inputdir/step1.root 
  if [ -f $PWD/step1.root ] ; then
      rm $PWD/step1.root
  fi
  ln -s $inputdir/step1_${i}.root $PWD/step1.root 
  #cmsRun $pydir/CellHitSum_cfg.py #-n 4
  cmsRun $pydir/SimHit_cfg.py #-n 4
  mv geantoutput.root geantoutput_${i}.root
done

ls $PWD/geantoutput_*.root > /tmp/idas/fl.txt
source ~/scripts/addhisto_file.sh /tmp/idas/fl.txt
#rm geantoutput_*.root
mv histo_merged.root geantoutput.root
