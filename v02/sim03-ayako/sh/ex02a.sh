#!/bin/sh
cd ..
unbuffer ./sim ICB/ex02/ 50 0.003 0.003 0.003 0 > log/ex02a_003.log
unbuffer ./sim ICB/ex02/ 50 0.004 0.004 0.004 0 > log/ex02a_004.log
unbuffer ./sim ICB/ex02/ 50 0.005 0.005 0.005 0 > log/ex02a_005.log
unbuffer ./sim ICB/ex02/ 50 0.006 0.006 0.006 0 > log/ex02a_006.log
unbuffer ./sim ICB/ex02/ 50 0.007 0.007 0.007 0 > log/ex02a_007.log
unbuffer ./sim ICB/ex02/ 50 0.008 0.008 0.008 0 > log/ex02a_008.log
unbuffer ./sim ICB/ex02/ 50 0.009 0.009 0.009 0 > log/ex02a_009.log
