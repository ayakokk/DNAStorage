#!/bin/sh
cd ../../
unbuffer ./sim ICB/ex04/ CSS/k03.txt 38 0.00 0.00 0.00 70 5 0 > log/D/e4k3.00.log
unbuffer ./sim ICB/ex04/ CSS/k03.txt 38 0.01 0.01 0.01 70 5 0 > log/D/e4k3.01.log
unbuffer ./sim ICB/ex04/ CSS/k03.txt 38 0.02 0.02 0.02 70 5 0 > log/D/e4k3.02.log
unbuffer ./sim ICB/ex04/ CSS/k03.txt 38 0.03 0.03 0.03 70 5 0 > log/D/e4k3.03.log
unbuffer ./sim ICB/ex04/ CSS/k03.txt 38 0.04 0.04 0.04 70 5 0 > log/D/e4k3.04.log
unbuffer ./sim ICB/ex04/ CSS/k03.txt 38 0.05 0.05 0.05 70 5 0 > log/D/e4k3.05.log
