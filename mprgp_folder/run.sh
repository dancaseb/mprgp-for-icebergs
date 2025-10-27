#!/bin/bash


ElmerSolver test-cases/case.sif

cd mprgp

# Step 1: Extract data from ParaView
pvpython ./save_paraview_data.py ../square1/case_default_solver_t0001.vtu output_default_solver.csv

# Step 2: Run Octave script that uses the CSV data
octave ./solve_LimitTemperature.m
