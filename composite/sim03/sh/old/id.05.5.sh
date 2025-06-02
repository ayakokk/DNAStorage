#!/bin/sh
cd ..
unbuffer ./sim ICB/ex01/ 50 0.01 0.01 0.00 5 5 0 > log/id010.05.5.000
unbuffer ./sim ICB/ex01/ 50 0.02 0.02 0.00 5 5 0 > log/id020.05.5.000
unbuffer ./sim ICB/ex01/ 50 0.03 0.03 0.00 5 5 0 > log/id030.05.5.000
unbuffer ./sim ICB/ex01/ 50 0.04 0.04 0.00 5 5 0 > log/id040.05.5.000
unbuffer ./sim ICB/ex01/ 50 0.05 0.05 0.00 5 5 0 > log/id050.05.5.000
