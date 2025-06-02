#!/bin/sh
cd ..
unbuffer ./SimInner 0.010 0.010 200 ../search03/ICB/G12.bin Gtable/4/010_010.bin 1 0 > log/G12/010.log
unbuffer ./SimInner 0.020 0.010 200 ../search03/ICB/G12.bin Gtable/4/020_010.bin 1 0 > log/G12/020.log
unbuffer ./SimInner 0.030 0.010 200 ../search03/ICB/G12.bin Gtable/4/030_010.bin 1 0 > log/G12/030.log
unbuffer ./SimInner 0.040 0.010 200 ../search03/ICB/G12.bin Gtable/4/040_010.bin 1 0 > log/G12/040.log
unbuffer ./SimInner 0.050 0.010 200 ../search03/ICB/G12.bin Gtable/4/050_010.bin 1 0 > log/G12/050.log
