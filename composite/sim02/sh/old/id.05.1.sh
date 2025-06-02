#!/bin/sh
cd ..
unbuffer ./sim ICB/ex01/ 50 0.01 0.01 0.00 5 1 0 > log/id010.05.1.000
unbuffer ./sim ICB/ex01/ 50 0.02 0.02 0.00 5 1 0 > log/id020.05.1.000
unbuffer ./sim ICB/ex01/ 50 0.03 0.03 0.00 5 1 0 > log/id030.05.1.000
unbuffer ./sim ICB/ex01/ 50 0.04 0.04 0.00 5 1 0 > log/id040.05.1.000
