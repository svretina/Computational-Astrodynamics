#!/bin/sh

rm symplectic
rm series
rm *0.txt
rm time.txt
# Compiling the two programs

icc -o symplectic 2_body_symplectic_6.c -lm
icc -o series 2_body_power_series.c -lm


# running the programs for many values
for e in 0.2 0.9
do
    for dt in 0.01 0.001
    do
	echo "e=" $e "dt=" $dt
	./symplectic $(echo $(./kep_2_cart.py $e) $dt) &
	
	./series $(echo $(./kep_2_cart.py $e) $dt 20) &
    done
done
