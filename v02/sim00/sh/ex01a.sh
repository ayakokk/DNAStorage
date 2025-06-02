#!/bin/sh
cd ..
unbuffer ./sim ICB/ex01/ 50 0.004 0.004 0.004 0 > log/ex01a_004.log
unbuffer ./sim ICB/ex01/ 50 0.005 0.005 0.005 0 > log/ex01a_005.log
unbuffer ./sim ICB/ex01/ 50 0.006 0.006 0.006 0 > log/ex01a_006.log
unbuffer ./sim ICB/ex01/ 50 0.007 0.007 0.007 0 > log/ex01a_007.log
unbuffer ./sim ICB/ex01/ 50 0.008 0.008 0.008 0 > log/ex01a_008.log
unbuffer ./sim ICB/ex01/ 50 0.009 0.009 0.009 0 > log/ex01a_009.log
