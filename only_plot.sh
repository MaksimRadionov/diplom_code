#!/bin/bash
#pkill gnuplot
#cat matrix_n_comp.txt matrix_comp.txt matrix_fcomp.txt> $1
#cat matrix_comp.txt > matrix.txt
#cat matrix_comp.txt > matrix.txt
file_name=$1
del_name=".txt"
jpgname=${file_name%$del_name}
echo $jpgname
gnuplot -p  << EOP
#set terminal wxt enhanced font 'Verdana,10' persist
set terminal wxt
#set terminal wxt 
#set term png
#set output "$jpgname.png"


set pm3d map
#set palette gray
#set palette model HSV
set samples 100; set isosamples 100


set title "$1"
#set cbrange [1035:1054]
splot "$1" matrix 
EOP
#eog data_1.jpg 
