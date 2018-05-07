#!/bin/bash
#pkill gnuplot_qt
rm ./matrix*.txt
make all_max && ./all_max < $1 $2 

gnuplot -p << EOP

#set terminal jpeg size 640,480
set terminal wxt
#set output "data_1.jpg"

splot "matrix_n_comp.txt" matrix with dots 
EOP

#eog data_1.jpg 
