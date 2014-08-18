#!/bin/bash

echo step0 zero beam spot size
./gen 100000 0 0 0 00 00 >& report.v3/step0.log

echo step1 nominal beam spot size and angular speread
./gen 100000 50 10 0.5 00 00 >& report.v3/step1.log

echo step1a nominal beam spot size zero angluar spread
./gen 100000 0 0 0.5 00 00 >& report.v3/step1a.log

echo step1b zero beam spot size and nominal angular spread
./gen 100000 50 10 0 00 00 >& report.v3/step1b.log

echo step2 nominal beam spot segmented pion tracker
./gen 100000 50 10 0.5 10 00 >& report.v3/step2.log

echo step2a nominal beam spot segmented pion tracker as well as strt detector
./gen 100000 50 10 0.5 11 00 >& report.v3/step2a.log

echo step3 nominal beam spot segmented pion tracker, MS on pion trackers 1 and 2 
./gen 100000 50 10 0.5 10 10 >& report.v3/step3.log

echo step3a nominal beam spot segmented pion tracker and start detectors, MS on pion trackers 1 and 2 
./gen 100000 50 10 0.5 11 10 >& report.v3/step3a.log

echo step3aI nominal beam spot segmented pion tracker and start detectors, start detector strip size halfed, MS on pion trackers 1 and 2 
./gen 100000 50 10 0.5 11 10 7 >& report.v3/step3aI.log

echo step3aII nominal beam spot segmented pion tracker and start detectors, start detector strip size quartered, MS on pion trackers 1 and 2 
./gen 100000 50 10 0.5 11 10 3.5 >& report.v3/step3aII.log

echo step3aIII nominal beam spot segmented pion tracker and start detectors, start detector strip size set to 1mm, MS on pion trackers 1 and 2 
./gen 100000 50 10 0.5 11 10 1 >& report.v3/step3aIII.log

echo step3b nominal beam spot segmented pion tracker, MS on pion trackers 1 and 2 and Start Detector
./gen 100000 50 10 0.5 10 11 >& report.v3/step3b.log

