#!/bin/sh
cd ../../
unbuffer ./sim ICB/ex04/ CSS/k03.txt 38 0.00 0.00 0.00 10 1 0 > log/A/e4k3.010.log
unbuffer ./sim ICB/ex04/ CSS/k03.txt 38 0.00 0.00 0.00 30 1 0 > log/A/e4k3.030.log
unbuffer ./sim ICB/ex04/ CSS/k03.txt 38 0.00 0.00 0.00 50 1 0 > log/A/e4k3.050.log
unbuffer ./sim ICB/ex04/ CSS/k03.txt 38 0.00 0.00 0.00 70 1 0 > log/A/e4k3.070.log
unbuffer ./sim ICB/ex04/ CSS/k03.txt 38 0.00 0.00 0.00 90 1 0 > log/A/e4k3.090.log
