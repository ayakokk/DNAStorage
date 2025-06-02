#!/bin/sh
cd ../../
unbuffer ./sim ICB/ex01/ CSS/k07.txt 50 0.00 0.00 0.00 20 1 0 > log/A/e1k7.020.log
unbuffer ./sim ICB/ex01/ CSS/k07.txt 50 0.00 0.00 0.00 40 1 0 > log/A/e1k7.040.log
unbuffer ./sim ICB/ex01/ CSS/k07.txt 50 0.00 0.00 0.00 60 1 0 > log/A/e1k7.060.log
unbuffer ./sim ICB/ex01/ CSS/k07.txt 50 0.00 0.00 0.00 80 1 0 > log/A/e1k7.080.log
