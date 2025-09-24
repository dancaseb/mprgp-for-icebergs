#!/bin/bash


ElmerSolver case.sif

cd mprgp
# Step 1: Extract data from ParaView
pvpython ./save_paraview_data.py

# Step 2: Run Octave script that uses the CSV data
octave ./solve_LimitTemperature.m
