#!/bin/bash
echo "Launching Tests for Alberto Module"

python Check_Models.py

python Check_MultiThreading.py 1

python Check_MultiThreading.py 4

python EBCC_closed_loop.py 1

mv OutputFile.dat Output_1.dat
mv CR.dat CR_1.dat

python EBCC_closed_loop.py 4

mv OutputFile.dat Output_4.dat
mv CR.dat CR_4.dat
