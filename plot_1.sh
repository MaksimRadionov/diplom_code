#!/bin/bash
echo "Usage ./.plot_1.sh <file_name> <height> <wigth> <time>"
rm ./graph.txt
#pkill gnuplot
make script && ./script  < $1 $2 $3 $4
title="$1height${2}wigth$3"
mv ./graph.txt $1
gnuplot -p << EOP

#set terminal jpeg size 640,480
#set output "data_1.jpg"
set terminal wxt
set title "$title"
plot "graph2.txt" with lines,\
        "$1" with lines 
EOP
#eog data_1.jpg 
