#!/bin/sh
cd ..
unbuffer ./sim ICB/ex02/ 50 0.002 0.002 0.000 0 > log/ex02b_002.log
unbuffer ./sim ICB/ex02/ 50 0.003 0.003 0.000 0 > log/ex02b_003.log
unbuffer ./sim ICB/ex02/ 50 0.004 0.004 0.000 0 > log/ex02b_004.log
unbuffer ./sim ICB/ex02/ 50 0.005 0.005 0.000 0 > log/ex02b_005.log
unbuffer ./sim ICB/ex02/ 50 0.006 0.006 0.000 0 > log/ex02b_006.log
unbuffer ./sim ICB/ex02/ 50 0.007 0.007 0.000 0 > log/ex02b_007.log
unbuffer ./sim ICB/ex02/ 50 0.008 0.008 0.000 0 > log/ex02b_008.log
unbuffer ./sim ICB/ex02/ 50 0.009 0.009 0.000 0 > log/ex02b_009.log
