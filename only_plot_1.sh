#!/bin/bash
pkill gnuplot_qt
title="height${1}wigth$2"
gnuplot -p << EOP

#set terminal jpeg size 640,480
#set output "data_1.jpg"

#set label $title
set terminal wxt
plot "graph2.txt" with lines,\
        "graph.txt" with lines 
EOP
#eog data_1.jpg 
