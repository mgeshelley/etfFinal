FC	= gfortran
FFLAGS	= -Wall -O0
# LDFLAGS	= -llapack -lblas
EXEC	= pairing_gap
SRC	= pairing_gap.f90 routines.f90
OBJ	= $(SRC:.f90=.o)

all: $(EXEC)

$(EXEC): $(OBJ)
	$(FC) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) -o $@ -c $< $(FFLAGS)

clean: 
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC)