# Makefile for coeff2D

FC = gfortran
LD = gfortran
SRC = src
BIN = bin
F90FLAGS = -fbounds-check -Wall -fbacktrace -g
FFLAGS = -O3
LDFLAGS= -llapack -lblas
_OBJECTS = type_defs.o quad_element.o legendre_module.o problemsetup.o weights.o main.o
OBJECTS = $(patsubst %,$(BIN)/%,$(_OBJECTS))
EXECUTABLE = main.x
# EXECUTABLE = test_pis.x

.PHONY: clean
$(EXECUTABLE): $(OBJECTS)
	       $(LD) -o $(EXECUTABLE) $(OBJECTS) -I$(BIN)/ $(LDFLAGS)
$(BIN)/%.o : $(SRC)/%.f90
	   $(FC) $(F90FLAGS) -c -o $@ $< -J$(BIN)/
$(BIN)/%.o : $(SRC)/%.f
	   $(FC) $(FFLAGS) -c -o $@ $<
clean:
	rm -f $(OBJECTS) $(EXECUTABLE) $(BIN)/*.mod
