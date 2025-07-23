#!/bin/sh
cd ..
unbuffer ./sim ICB/ex02/ 50 0.000 0.002 0.000 0 > log/ex02c_002.log
unbuffer ./sim ICB/ex02/ 50 0.000 0.003 0.000 0 > log/ex02c_003.log
unbuffer ./sim ICB/ex02/ 50 0.000 0.004 0.000 0 > log/ex02c_004.log
unbuffer ./sim ICB/ex02/ 50 0.000 0.005 0.000 0 > log/ex02c_005.log
unbuffer ./sim ICB/ex02/ 50 0.000 0.006 0.000 0 > log/ex02c_006.log
unbuffer ./sim ICB/ex02/ 50 0.000 0.007 0.000 0 > log/ex02c_007.log
unbuffer ./sim ICB/ex02/ 50 0.000 0.008 0.000 0 > log/ex02c_008.log
unbuffer ./sim ICB/ex02/ 50 0.000 0.009 0.000 0 > log/ex02c_009.log
