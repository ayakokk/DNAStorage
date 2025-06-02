#!/bin/sh
cd ../../
unbuffer ./sim ICB/ex04/ CSS/k07.txt 38 0.00 0.00 0.00 10 1 0 > log/A/e4k7.010.log
unbuffer ./sim ICB/ex04/ CSS/k07.txt 38 0.00 0.00 0.00 30 1 0 > log/A/e4k7.030.log
unbuffer ./sim ICB/ex04/ CSS/k07.txt 38 0.00 0.00 0.00 50 1 0 > log/A/e4k7.050.log
unbuffer ./sim ICB/ex04/ CSS/k07.txt 38 0.00 0.00 0.00 70 1 0 > log/A/e4k7.070.log
unbuffer ./sim ICB/ex04/ CSS/k07.txt 38 0.00 0.00 0.00 90 1 0 > log/A/e4k7.090.log
