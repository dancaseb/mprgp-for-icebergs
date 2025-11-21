ALL: mprgp mprgp_mat_setup elmer

CC = mpicc
CFLAGS = -O2 -g -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PERMON_DIR)/include
LDFLAGS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(PERMON_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,$(PERMON_DIR)/$(PETSC_ARCH)/lib \
	-Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,$(PERMON_DIR)/$(PETSC_ARCH)/lib
LDLIBS  = -lpetsc -lpermon -lm

mprgp: mprgp.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS) $(LDLIBS)
mprgp_mat_setup: mprgp_mat_setup.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS) $(LDLIBS)
elmer: elmer.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS) $(LDLIBS)
clean:
	rm -f mprgp
	rm -f mprgp_mat_setup
	rm -f elmer