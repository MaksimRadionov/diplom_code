#!/bin/bash
rm ./graph.txt
pkill gnuplot_qt
make test && ./test  
title="height${1}wigth$2"
gnuplot -p << EOP

#set terminal jpeg size 640,480
#set output "data_1.jpg"

#set label $title
plot "graph2.txt" with lines,\
        "graph.txt" with points 
EOP
#eog data_1.jpg 
