#!/bin/sh
cd ../../
unbuffer ./sim ICB/ex04/ CSS/k07.txt 38 0.00 0.00 0.00 30 2 0 > log/B/30_2.00.log
unbuffer ./sim ICB/ex04/ CSS/k07.txt 38 0.01 0.01 0.01 30 2 0 > log/B/30_2.01.log
unbuffer ./sim ICB/ex04/ CSS/k07.txt 38 0.02 0.02 0.02 30 2 0 > log/B/30_2.02.log
unbuffer ./sim ICB/ex04/ CSS/k07.txt 38 0.03 0.03 0.03 30 2 0 > log/B/30_2.03.log
unbuffer ./sim ICB/ex04/ CSS/k07.txt 38 0.04 0.04 0.04 30 2 0 > log/B/30_2.04.log
unbuffer ./sim ICB/ex04/ CSS/k07.txt 38 0.05 0.05 0.05 30 2 0 > log/B/30_2.05.log
