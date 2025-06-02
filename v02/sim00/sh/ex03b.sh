#!/bin/sh
cd ..
unbuffer ./sim ICB/ex03/ 50 0.002 0.002 0.000 0 > log/ex03b_002.log
unbuffer ./sim ICB/ex03/ 50 0.003 0.003 0.000 0 > log/ex03b_003.log
unbuffer ./sim ICB/ex03/ 50 0.004 0.004 0.000 0 > log/ex03b_004.log
unbuffer ./sim ICB/ex03/ 50 0.005 0.005 0.000 0 > log/ex03b_005.log
unbuffer ./sim ICB/ex03/ 50 0.006 0.006 0.000 0 > log/ex03b_006.log
unbuffer ./sim ICB/ex03/ 50 0.007 0.007 0.000 0 > log/ex03b_007.log
unbuffer ./sim ICB/ex03/ 50 0.008 0.008 0.000 0 > log/ex03b_008.log
unbuffer ./sim ICB/ex03/ 50 0.009 0.009 0.000 0 > log/ex03b_009.log
