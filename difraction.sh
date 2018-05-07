#!/bin/bash
rm init_signal.txt difr_signal.txt
make all_max && ./all_max
pkill gnuplot_qt
gnuplot -p << EOP

set terminal wxt
plot  "init_signal.txt" with lines,\
        "difr_signal.txt" with lines 
EOP

