# So that bash commands involving tee are read correctly by make
SHELL:=/bin/bash

# Add .f90, .o to list of acceptable suffixes
.SUFFIXES: .f90 .o $(SUFFIXES)

FC	= gfortran
# Use 2008 Fortran standard, and allow "dfloat" function (GNU extension)
FFLAGS	= -std=f2008 -fall-intrinsics
# All warnings and checks, showing maximum of one error
WFLAGS_1	= -Wall -Wextra -Wconversion -fcheck=all -fmax-errors=1
# Remove annoying warnings in "boundary.f90", coming from unused variables, from unused labels, and from real comparisons
WFLAGS_2	= -Wall -Wextra -Wconversion -Wno-unused-variable -Wno-unused-label -Wno-compare-reals -fcheck=all -fmax-errors=1
# LDFLAGS	= -llapack -lblas

EXEC	= etf
EXEC_OBJ= $(EXEC).o

# .f90 files stored in SRCDIR
SRCDIR	= srcs
# Two separate lists of source files, for compilation with different warning flags
# wildcard selects all .f90 files
# filter-out removes boundary.f90
SRC_1	:= $(filter-out $(SRCDIR)/boundary.f90, $(wildcard $(SRCDIR)/*.f90))
# addprefix adds SRCDIR prefix
SRC_2	:= $(addprefix $(SRCDIR)/,boundary.f90)
SRC	:= $(SRC_1) $(SRC_2)

# .o files stored in MODDIR
MODDIR	= mods
# Two separate lists of object files
# notdir remove SRCDIR prefix
# : and = is for replacing .f90 with .o suffixes
# addprefix adds MODDIR prefix
OBJ_1	:= $(addprefix $(MODDIR)/,$(notdir $(SRC_1:.f90=.o)))
OBJ_2	:= $(addprefix $(MODDIR)/,$(notdir $(SRC_2:.f90=.o)))
OBJ	:= $(OBJ_1) $(OBJ_2)


all: $(EXEC)

# Assign warnings to "WFLAGS", if needed
$(OBJ_1): WFLAGS := $(WFLAGS_1)
$(OBJ_2): WFLAGS := $(WFLAGS_2)

# Add warning flags if "make debug"
debug: FFLAGS += $(WFLAGS)
debug: clean $(EXEC)

# Add -O3 for speedup if "make fast"
fast: FFLAGS += -O3
fast: clean $(EXEC)

# Add -fopenmp for parallel if "make fast_parallel"
fast_parallel: FFLAGS += -fopenmp
fast_parallel: fast $(EXEC)

# Add -pg flag for profiling
profile: FFLAGS = -pg
profile: clean $(EXEC)

# -p means ignore warning if directory already exists
$(SRCDIR):
	mkdir -p $@

$(MODDIR):
	mkdir -p $@

# All object files required for EXEC, use $^
$(EXEC): $(OBJ)
	$(FC) -o $@ $^ $(FFLAGS)

# $(MODDIR)/%.o means .o files created are stored in MODDIR
# % with $< means each file in %.f90 will be compiled separately
# -c means no linking of .o files occurs
# -J$| stores .mod files in MODDIR
$(MODDIR)/%.o: $(SRCDIR)/%.f90 | $(MODDIR)
	$(FC) -o $@ -c $< -J$| $(FFLAGS)

clean: 
	rm -rf $(MODDIR)/*.o
	rm -rf $(MODDIR)/*.mod

mrproper: clean
	rm -rf $(EXEC)

# Print all variables for debugging
debug_make:
	@echo "FC = $(FC)"
	@echo "FFLAGS = $(FFLAGS)"
	@echo "WFLAGS_1 = $(WFLAGS_1)"
	@echo "WFLAGS_2 = $(WFLAGS_2)"
	@echo "MODDIR = $(MODDIR)"
	@echo "SRCDIR = $(SRCDIR)"
	@echo "EXEC = $(EXEC)"
	@echo "EXEC_OBJ = $(EXEC_OBJ)"
	@echo "SRC_1 = $(SRC_1)"
	@echo "SRC_2 = $(SRC_2)"
	@echo "SRC = $(SRC)"
	@echo "OBJ_1 = $(OBJ_1)"
	@echo "OBJ_2 = $(OBJ_2)"
	@echo "OBJ = $(OBJ)"
	@echo "test = $(MODDIR)/%.o"

# Make, then run, copying all screen output to "output.out"
# @ at beginning means don't print any commands that are run
${EXEC:%=run_%}: run_%: $(EXEC)
	@./$< |& tee output.out

# Run with no screen output
run_silent: $(EXEC)
	@./$< > output.out


.PHONY: all debug clean mrproper debug_make run_%


# Dependencies
EXEC_OBJ : $(OBJ) # Main program depends on all modules
$(MODDIR)/routines.o: $(MODDIR)/parameters.o
$(MODDIR)/boundary.o: $(MODDIR)/parameters.o
$(MODDIR)/strutinsky.o: $(MODDIR)/parameters.o $(MODDIR)/routines.o $(MODDIR)/boundary.o 
$(MODDIR)/pairing.o: $(MODDIR)/parameters.o $(MODDIR)/routines.o $(MODDIR)/strutinsky.o $(MODDIR)/test.o
$(MODDIR)/test.o: $(MODDIR)/parameters.o $(MODDIR)/routines.o $(MODDIR)/boundary.o $(MODDIR)/strutinsky.o