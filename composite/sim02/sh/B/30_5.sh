#!/bin/sh
cd ../../
unbuffer ./sim ICB/ex04/ CSS/k07.txt 38 0.00 0.00 0.00 30 5 0 > log/B/30_5.00.log
unbuffer ./sim ICB/ex04/ CSS/k07.txt 38 0.01 0.01 0.01 30 5 0 > log/B/30_5.01.log
unbuffer ./sim ICB/ex04/ CSS/k07.txt 38 0.02 0.02 0.02 30 5 0 > log/B/30_5.02.log
unbuffer ./sim ICB/ex04/ CSS/k07.txt 38 0.03 0.03 0.03 30 5 0 > log/B/30_5.03.log
unbuffer ./sim ICB/ex04/ CSS/k07.txt 38 0.04 0.04 0.04 30 5 0 > log/B/30_5.04.log
unbuffer ./sim ICB/ex04/ CSS/k07.txt 38 0.05 0.05 0.05 30 5 0 > log/B/30_5.05.log
