#!/bin/sh
cd ..
unbuffer ./SimInner 0.020 0.010 200 ../search03/ICB/A12.bin Gtable/5/020_010.bin 1 0 > log/A12/020.log
unbuffer ./SimInner 0.030 0.010 200 ../search03/ICB/A12.bin Gtable/5/030_010.bin 1 0 > log/A12/030.log
unbuffer ./SimInner 0.040 0.010 200 ../search03/ICB/A12.bin Gtable/5/040_010.bin 1 0 > log/A12/040.log
unbuffer ./SimInner 0.050 0.010 200 ../search03/ICB/A12.bin Gtable/5/050_010.bin 1 0 > log/A12/050.log
