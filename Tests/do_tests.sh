#!/bin/bash
echo "Launching Tests for Alberto Module"

python Check_Models.py

mpirun -np 4 python Check_MultiThreading.py

