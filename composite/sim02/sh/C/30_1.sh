#!/bin/sh
cd ../../
unbuffer ./sim ICB/ex04/ CSS/k07.txt 38 0.00 0.00 0.00 30 1 0 > log/C/30_1.00.log
unbuffer ./sim ICB/ex04/ CSS/k07.txt 38 0.01 0.01 0.00 30 1 0 > log/C/30_1.01.log
unbuffer ./sim ICB/ex04/ CSS/k07.txt 38 0.02 0.02 0.00 30 1 0 > log/C/30_1.02.log
unbuffer ./sim ICB/ex04/ CSS/k07.txt 38 0.03 0.03 0.00 30 1 0 > log/C/30_1.03.log
unbuffer ./sim ICB/ex04/ CSS/k07.txt 38 0.04 0.04 0.00 30 1 0 > log/C/30_1.04.log
unbuffer ./sim ICB/ex04/ CSS/k07.txt 38 0.05 0.05 0.00 30 1 0 > log/C/30_1.05.log
