TARGET = qutraj_run

F90 = gfortran
F90FLAGS =
LDFLAGS = -lblas
F77 = gfortran
F77FLAGS =

PYFFLAGS = --overwrite-signature
F2PYBUILDFLAGS = -DF2PY_REPORT_ON_ARRAY_COPY=1

ZVODE = zvode zcopy zgbfa zgbsl zgefa zgesl
ZVODE_OBJS = $(ZVODE:=.o)
ZVODE_SRCS = $(ZVODE:=.f)

MODS = qutraj_general
MODS += qutraj_hilbert
MODS += qutraj_solver

OBJS = $(MODS:=.o)
SRCS += $(MODS:=.f90)

#MODS += $(ZVODE)
#OBJS += $(ZVODE_OBJS)
#SRCS += $(ZVODE_SRCS)

.PHONY: all clean

all:
	@echo Usage:
	@echo make pyf - makes .pyf file, do this first.
	@echo make f2pybuild - compiles the module and all dependencies.

pyf:
ifneq "$(wildcard qutraj_run.pyf )" ""
	-cp $(TARGET:=.pyf) $(TARGET:=.pyf.backup)
endif
	f2py -m $(TARGET) -h $(TARGET:=.pyf) $(PYFFLAGS) $(TARGET:=.f90)

#
# need to manually edit .pyf file to remove private subroutines,
# and give default values for optional arguments
#

f2pybuild: $(OBJS)
	f2py $(F2PYBUILDFLAGS) -c $(TARGET:=.pyf) $(TARGET:=.f90) $(SRCS)

$(OBJS): $(ZVODE_OBJS)

#f90mcsolve: $(OBJS) f90mcsolve.o
	#$(F90) $(F90FLAGS) -c $(OBJS) $(@:=.o)

clean:
	$(RM) *.o
	$(RM) *.mod
	$(RM) *.so


.SUFFIXES: $(SUFFIXES) .f90 .f

# Need to override standard rule for mod files since make thinks they are modula-2 files
%.o: %.mod
.f90.o:
# Rule on how to convert from .f90 to .o
# $< macro refers to the name of the related file that caused the action.
	$(F90) $(F90FLAGS) -c $< 
.f90.mod:
	$(F90) $(F90FLAGS) -c $< 
.f.o:
	$(F77) $(F77FLAGS) -c $< 

#
# Dependencies
#

qutraj_run.o: qutraj_run.f90 qutraj_solver.o qutraj_general.o
qutraj_solver.o: qutraj_solver.f90 qutraj_hilbert.o qutraj_general.o
qutraj_hilbert.o: qutraj_hilbert.f90 qutraj_general.o
qutraj_general.o: qutraj_general.f90

