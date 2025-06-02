#!/bin/sh
cd ../../
unbuffer ./sim ICB/ex01/ CSS/k03.txt 50 0.00 0.00 0.00 30 1 0 > log/A/e1k3.030.log
unbuffer ./sim ICB/ex01/ CSS/k03.txt 50 0.00 0.00 0.00 50 1 0 > log/A/e1k3.050.log
unbuffer ./sim ICB/ex01/ CSS/k03.txt 50 0.00 0.00 0.00 70 1 0 > log/A/e1k3.070.log
unbuffer ./sim ICB/ex01/ CSS/k03.txt 50 0.00 0.00 0.00 90 1 0 > log/A/e1k3.090.log
