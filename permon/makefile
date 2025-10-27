ALL: mprgp

CC = mpicc
CFLAGS = -O2 -g -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PERMON_DIR)/include
LDFLAGS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(PERMON_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,$(PERMON_DIR)/$(PETSC_ARCH)/lib \
	-Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,$(PERMON_DIR)/$(PETSC_ARCH)/lib
LDLIBS  = -lpetsc -lpermon -lm

mprgp: mprgp.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS) $(LDLIBS)
clean:
	rm -f mprgp
