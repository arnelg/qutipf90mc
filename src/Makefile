TARGET = f90mcsolve
#TARGET = f2pytest

F90 = gfortran
F90FLAGS = -O3
LDFLAGS =

PYFFLAGS = --overwrite-signature
F2PYBUILDFLAGS = -DF2PY_REPORT_ON_ARRAY_COPY=1

MODS = solver

OBJS = $(MODS:=.o)
SRCS += $(MODS:=.f90)

.PHONY: all clean

all:
	@echo Usage:
	@echo make pyf - makes .pyf file, do this first.
	@echo make f2pybuild - compiles the module and all dependencies.

pyf:
	f2py -m $(TARGET) -h $(TARGET:=.pyf) $(PYFFLAGS) $(TARGET:=.f90)
	
#
# need to manually edit .pyf file to remove private subroutines.
#

f2pybuild: $(OBJS)
	f2py $(F2PYBUILDFLAGS) -c $(TARGET:=.pyf) $(TARGET:=.f90) $(SRCS)

#f90mcsolve: $(OBJS) f90mcsolve.o
	#$(F90) $(F90FLAGS) -c $(OBJS) $(@:=.o)

clean:
	$(RM) *.o
	$(RM) *.mod
	$(RM) *.so


.SUFFIXES: $(SUFFIXES) .f90

# Need to override standard rule for mod files since make thinks they are modula-2 files
%.o: %.mod
.f90.o:
# Rule on how to convert from .f90 to .o
# $< macro refers to the name of the related file that caused the action.
	$(F90) $(F90FLAGS) -c $< 
.f90.mod:
	$(F90) $(F90FLAGS) -c $< 

#
# Dependencies
#

f90mcsolve.o: f90mcsolve.f90 solver.o 
solver.o: solver.f90
