#!/bin/bash


ElmerSolver case.sif
ElmerSolver case_my_solver.sif

cd mprgp

# Step 1: Extract data from ParaView
pvpython ./save_paraview_data.py  "../square/case_default_solver_t0001.vtu" "output_default_solver.csv"
pvpython ./save_paraview_data.py  "../square/case_my_solver_t0001.vtu" "output_my_solver.csv"


octave ./compare_x.m



