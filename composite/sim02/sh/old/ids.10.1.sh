#!/bin/sh
cd ..
unbuffer ./sim ICB/ex01/ 50 0.01 0.01 0.01 10 1 0 > log/ids010.10.1.000
unbuffer ./sim ICB/ex01/ 50 0.02 0.02 0.02 10 1 0 > log/ids020.10.1.000
unbuffer ./sim ICB/ex01/ 50 0.03 0.03 0.03 10 1 0 > log/ids030.10.1.000
unbuffer ./sim ICB/ex01/ 50 0.04 0.04 0.04 10 1 0 > log/ids040.10.1.000
