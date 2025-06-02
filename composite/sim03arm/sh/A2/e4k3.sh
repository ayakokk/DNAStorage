#!/bin/sh
cd ../../
unbuffer ./sim ICB/ex04/ CSS/k03.txt 38 0.00 0.00 0.00 20 1 0 > log/A/e4k3.020.log
unbuffer ./sim ICB/ex04/ CSS/k03.txt 38 0.00 0.00 0.00 40 1 0 > log/A/e4k3.040.log
unbuffer ./sim ICB/ex04/ CSS/k03.txt 38 0.00 0.00 0.00 60 1 0 > log/A/e4k3.060.log
unbuffer ./sim ICB/ex04/ CSS/k03.txt 38 0.00 0.00 0.00 80 1 0 > log/A/e4k3.080.log
