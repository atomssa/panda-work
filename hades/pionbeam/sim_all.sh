#!/bin/bash

nevt=10000
if [ "$1" == 1 ];
then
    echo arg = 1
    echo step0 zero beam spot size
    ./gen ${nevt} 0 0 0 00 00 >& output/report.v4/step0.log
    echo step1 nominal beam spot size and angular speread
    ./gen ${nevt} 50 10 0.5 00 00 >& output/report.v4/step1.log
    echo step1a nominal beam spot size zero angluar spread
    ./gen ${nevt} 0 0 0.5 00 00 >& output/report.v4/step1a.log
    echo step1b zero beam spot size and nominal angular spread
    ./gen ${nevt} 50 10 0 00 00 >& output/report.v4/step1b.log
fi

if [ "$1" == 2 ];
then
    echo arg = 2
    # seg on poin tracker 1 and 2  only
    echo step2 nominal beam spot segmented pion tracker no MS
    ./gen ${nevt} 50 10 0.5 10 00 >& output/report.v4/step2.log

    echo step2a nominal beam spot segmented pion tracker, MS on pion trackers 1
    ./gen ${nevt} 50 10 0.5 10 10 >& output/report.v4/step2a.log

    echo step2b nominal beam spot segmented pion tracker, MS on pion trackers 1 and 2
    ./gen ${nevt} 50 10 0.5 10 11 >& output/report.v4/step2b.log
fi

if [ "$1" == 3 ];
then
    echo arg = 3
    # seg on poin trackers and start detector
    echo step3 nominal beam spot segmented pion tracker 1 and 2 as well as strt detector
    ./gen ${nevt} 50 10 0.5 11 00 >& output/report.v4/step3.log

    echo step3a nominal beam spot segmented pion tracker as well as strt detector, MS on pion trackers 1
    ./gen ${nevt} 50 10 0.5 11 10 >& output/report.v4/step3a.log

    echo step3b nominal beam spot segmented pion tracker as well as strt detector, MS on pion trackers 1 and 2
    ./gen ${nevt} 50 10 0.5 11 11 >& output/report.v4/step3b.log
fi

#echo step3aI nominal beam spot segmented pion tracker and start detectors, start detector strip size halfed, MS on pion trackers 1 and 2
#./gen 100000 50 10 0.5 11 10 7 >& output/report.v4/step3aI.log
#
#echo step3aII nominal beam spot segmented pion tracker and start detectors, start detector strip size quartered, MS on pion trackers 1 and 2
#./gen 100000 50 10 0.5 11 10 3.5 >& output/report.v4/step3aII.log
#
#echo step3aIII nominal beam spot segmented pion tracker and start detectors, start detector strip size set to 1mm, MS on pion trackers 1 and 2
#./gen 100000 50 10 0.5 11 10 1 >& output/report.v4/step3aIII.log

#echo step3b nominal beam spot segmented pion tracker, MS on pion trackers 1 and 2 and Start Detector
#./gen 100000 50 10 0.5 11 11 >& output/report.v4/step3b.log
#
#echo step3c nominal beam spot segmented pion tracker and start detector, MS on pion trackers 1 and 2 and Start Detector
#./gen 100000 50 10 0.5 11 11 >& output/report.v4/step3c.log
