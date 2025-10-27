# Make sure to call the ElmerSolver to get your data files first (for testing)
# call make

# Convert elmer data files to PETSc binary format
octave --quiet mat_elmer_to_petsc.m linsys_a.dat linsys_a.bin
octave --quiet vec_elmer_to_petsc.m linsys_b.dat linsys_b.bin
octave --quiet vec_elmer_to_petsc.m c.dat c.bin

# Run the MPRGP solver
./mprgp -a linsys_a.bin -b linsys_b.bin -c c.bin

octave --quiet vec_petsc_to_elmer.m x_output.bin x.dat

octave --quiet compare_permon_elmer.m

