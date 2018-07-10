#!/bin/bash
echo "Launching Tests for Alberto Module"

python Check_Models.py

mpirun -np 1 python Check_MultiThreading.py 1

mpirun -np 4 python Check_MultiThreading.py 4

python Remove_Empty.py

diff -u  Weights_1-* Weights_4-* > Diff_Weights.csv

if [[ -s Diff_Weights.csv ]]; then echo "WARNING! Weight Files are different!"; else echo "OK"; fi

mpirun -np 1 python EBCC_closed_loop.py 1

mv OutputFile.dat Output_1.dat
mv CR.dat CR_1.dat

mpirun -np 4 python EBCC_closed_loop.py 4

mv OutputFile.dat Output_4.dat
mv CR.dat CR_4.dat
