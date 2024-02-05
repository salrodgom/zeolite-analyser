# Salvador R. G. Balestra, 2020
# Makefile for main.f90
FFLAGS = -ffree-form -march=native -O2 -fbounds-check
PROGRAM = zeolite_analyser.exe
install:
	gfortran main.f90 -o $(PROGRAM) $(FFLAGS)
all:
	gfortran main.f90 -o $(PROGRAM) $(FFLAGS)
	./$(PROGRAM)
	rm -rf $(PROGRAM) *.o *.mod
clean:
	rm -rf $(PROGRAM) *.o *.mod
