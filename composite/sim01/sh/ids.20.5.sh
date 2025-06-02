#!/bin/sh
cd ..
unbuffer ./sim ICB/ex01/ 50 0.01 0.01 0.01 20 5 0 > log/ids010.20.5.000
unbuffer ./sim ICB/ex01/ 50 0.02 0.02 0.02 20 5 0 > log/ids020.20.5.000
unbuffer ./sim ICB/ex01/ 50 0.03 0.03 0.03 20 5 0 > log/ids030.20.5.000
unbuffer ./sim ICB/ex01/ 50 0.04 0.04 0.04 20 5 0 > log/ids040.20.5.000
unbuffer ./sim ICB/ex01/ 50 0.05 0.05 0.05 20 5 0 > log/ids050.20.5.000
