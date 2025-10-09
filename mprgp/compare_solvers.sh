#!/bin/bash

# compile the custom solver
elmerf90 -o MyPDE2.so MyPDE2.F90
if [ $? -eq 0 ]; then
    echo "Compilation successful!"
else
    echo "Compilation failed!"
    exit 1
fi

# Step 1: Run ElmerSolver with different solvers
rm -f output_default_solver.csv
rm -f output_my_solver.csv
rm -f ../square/case_default_solver_t0001.vtu
rm -f ../square/case_my_solver_t0001.vtu

ElmerSolver test-cases/case.sif
ElmerSolver test-cases/case_my_solver.sif

cd mprgp

# Step 2: Extract data from ParaView
pvpython ./save_paraview_data.py  "../square/case_default_solver_t0001.vtu" "output_default_solver.csv"
pvpython ./save_paraview_data.py  "../square/case_my_solver_t0001.vtu" "output_my_solver.csv"


# Step 3: Compare results using Octave
octave ./compare_x.m



