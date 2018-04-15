.SUFFIXES: .f .o
.SUFFIXES: .f90 .o
.SUFFIXES: .f90 .mod

BINDIR = ./
FC = gfortran

FFLAGS = -O1  -fbacktrace -g -Wall -fbounds-check
FINC   = 
LDFLAGS= -llapack -lblas

LD = $(FC) 

# Object files 


OBJS_LDG_B_2D = ldg_burger_2d.o type_defs.o legfuns.o problem_burger_2d.o lglnodes.o 

#
# 
.PHONY:copy

$(BINDIR)/ldg_burger_2d.x: $(OBJS_LDG_B_2D)
	$(LD) $(LDFLAGS) -o $@ $(OBJS_LDG_B_2D) 


.f.o:	
	$(FC) -c $(FFLAGS) $(FINC) $<
.f90.o:	
	$(FC) -c $(FFLAGS) $(FINC) $<
.f90.mod:	
	$(FC) -c $(FFLAGS) $(FINC) $<

# The dependencies between modules
ldg_burger_2d.f90: type_defs.mod problem_burger_2d.mod legfuns.mod lglnodes.o	
legfuns.f90: type_defs.mod
leglnodes.f90: type_defs.mod
problem_burger_2d.f90: type_defs.mod

clean:
	rm -f *.o *.a *.mod *.txt
distclean:
	rm -f *.o *.a *.mod $(BASE)
