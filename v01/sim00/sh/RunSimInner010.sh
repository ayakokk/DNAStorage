#!/bin/sh
cd ..
unbuffer ./SimInner 0.010 0.010 200 ICB/1/5_08_2.bin Gtable/5/010_010.bin 0 > tmp/010_010.log
unbuffer ./SimInner 0.015 0.010 200 ICB/1/5_08_2.bin Gtable/5/015_010.bin 0 > tmp/015_010.log
unbuffer ./SimInner 0.020 0.010 200 ICB/1/5_08_2.bin Gtable/5/020_010.bin 0 > tmp/020_010.log
unbuffer ./SimInner 0.025 0.010 200 ICB/1/5_08_2.bin Gtable/5/025_010.bin 0 > tmp/025_010.log
unbuffer ./SimInner 0.030 0.010 200 ICB/1/5_08_2.bin Gtable/5/030_010.bin 0 > tmp/030_010.log
unbuffer ./SimInner 0.035 0.010 200 ICB/1/5_08_2.bin Gtable/5/035_010.bin 0 > tmp/035_010.log
unbuffer ./SimInner 0.040 0.010 200 ICB/1/5_08_2.bin Gtable/5/040_010.bin 0 > tmp/040_010.log
unbuffer ./SimInner 0.045 0.010 200 ICB/1/5_08_2.bin Gtable/5/045_010.bin 0 > tmp/045_010.log
unbuffer ./SimInner 0.050 0.010 200 ICB/1/5_08_2.bin Gtable/5/050_010.bin 0 > tmp/050_010.log
