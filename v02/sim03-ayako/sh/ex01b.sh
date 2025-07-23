#!/bin/sh
cd ..
unbuffer ./sim ICB/ex01/ 50 0.002 0.002 0.000 0 > log/ex01b_002.log
unbuffer ./sim ICB/ex01/ 50 0.003 0.003 0.000 0 > log/ex01b_003.log
unbuffer ./sim ICB/ex01/ 50 0.004 0.004 0.000 0 > log/ex01b_004.log
unbuffer ./sim ICB/ex01/ 50 0.005 0.005 0.000 0 > log/ex01b_005.log
unbuffer ./sim ICB/ex01/ 50 0.006 0.006 0.000 0 > log/ex01b_006.log
unbuffer ./sim ICB/ex01/ 50 0.007 0.007 0.000 0 > log/ex01b_007.log
unbuffer ./sim ICB/ex01/ 50 0.008 0.008 0.000 0 > log/ex01b_008.log
unbuffer ./sim ICB/ex01/ 50 0.009 0.009 0.000 0 > log/ex01b_009.log
