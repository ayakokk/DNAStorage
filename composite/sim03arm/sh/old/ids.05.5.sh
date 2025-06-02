#!/bin/sh
cd ..
unbuffer ./sim ICB/ex01/ 50 0.01 0.01 0.01 5 5 0 > log/ids010.05.5.000
unbuffer ./sim ICB/ex01/ 50 0.02 0.02 0.02 5 5 0 > log/ids020.05.5.000
unbuffer ./sim ICB/ex01/ 50 0.03 0.03 0.03 5 5 0 > log/ids030.05.5.000
unbuffer ./sim ICB/ex01/ 50 0.04 0.04 0.04 5 5 0 > log/ids040.05.5.000
