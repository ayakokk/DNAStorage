#!/bin/sh
cd ..
unbuffer ./sim ICB/ex04/ 38 0.002 0.002 0.002 0 > log/ex04a_002.log
unbuffer ./sim ICB/ex04/ 38 0.003 0.003 0.003 0 > log/ex04a_003.log
unbuffer ./sim ICB/ex04/ 38 0.004 0.004 0.004 0 > log/ex04a_004.log
unbuffer ./sim ICB/ex04/ 38 0.005 0.005 0.005 0 > log/ex04a_005.log
unbuffer ./sim ICB/ex04/ 38 0.006 0.006 0.006 0 > log/ex04a_006.log
unbuffer ./sim ICB/ex04/ 38 0.007 0.007 0.007 0 > log/ex04a_007.log
unbuffer ./sim ICB/ex04/ 38 0.008 0.008 0.008 0 > log/ex04a_008.log
unbuffer ./sim ICB/ex04/ 38 0.009 0.009 0.009 0 > log/ex04a_009.log
