# Universal Makefile for sceptic3D


# Set shell to bash (default is sh)
SHELL := /bin/bash

# Set compilers
G77 := $(shell ./setcomp f77)
G77nonmpi := $(shell ./setcomp f77 nonmpi)
G90 := $(shell ./setcomp f90)
G90nonmpi := $(shell ./setcomp f90 nonmpi)

# Set Xlib location
DIRXLIB := $(shell ./setxlib)
# Note that this may not work properly, and is not needed
#   if X11 and Xt are in /usr/lib(64)

# Accis lib location unless already set
DIRACCIS ?= ./accis


# Libraries and options to pass to linker
LIB := -L$(DIRXLIB) -L$(DIRACCIS) -laccisX -lXt -lX11 
# Show time and memory usage (debugging)
LIB += -Wl,-stats


# Options to pass to compiler
OPTCOMP := -I.
# Show all warnings exept unused variables
OPTCOMP += -Wall -Wno-unused-variable
# Enable optimization
OPTCOMP += -O2
# Save debugging info
OPTCOMP += -g
# Do bounds check (debugging)
#OPTCOMP += -ffortran-bounds-check
# Save profiling information (debugging)
OPTCOMP += -pg

# Options to pass to compiler for MPI version
OPTCOMPMPI := $(OPTCOMP)
# Enable MPI by defining 'MPI' for pre-compiler
OPTCOMPMPI += -DMPI


# Objects common to all versions of sceptic3D
OBJ := initiate.o \
       advancing.o \
       randc.o \
       randf.o \
       diags.o \
       outputs.o \
       chargefield.o \
       stringsnames.o \
       rhoinfcalc.o \
       shielding3D.o
# Reinjection related objects
OBJ += orbitinject.o \
       extint.o \
       maxreinject.o

# Objects for MPI version of sceptic3D
OBJMPI := $(OBJ) \
          cg3dmpi.o \
          mpibbdy.o \
          shielding3D_par.o

# Default target is serial sceptic3D without HDF support
sceptic3D : sceptic3D.F piccom.f $(OBJ) ./accis/libaccisX.a
	$(G77) $(OPTCOMP) -o sceptic3D sceptic3D.F $(OBJ) $(LIB)

# sceptic3D with MPI
sceptic3Dmpi : sceptic3D.F piccom.f piccomcg.f $(OBJMPI) ./accis/libaccisX.a
	$(G77) $(OPTCOMPMPI) -o sceptic3Dmpi sceptic3D.F $(OBJMPI) $(LIB)


# Other rules
./accis/libaccisX.a : ./accis/*.f
	make -C accis

orbitint : orbitint.f coulflux.o $(OBJ) ./accis/libaccisX.a
	$(G77) $(OPTCOMP) -o orbitint orbitint.f $(OBJ) coulflux.o $(LIB)

coulflux.o : tools/coulflux.f
	$(G77) -c $(OPTCOMP) tools/coulflux.f

fvinjecttest : fvinjecttest.F fvinject.o reinject.o initiate.o advancing.o chargefield.o randf.o fvcom.f
	$(G77)  -o fvinjecttest $(OPTCMOP) fvinjecttest.F fvinject.o reinject.o initiate.o advancing.o chargefield.o randf.o  $(LIB)

fvinject.o : fvinject.f fvcom.f piccom.f
	$(G77) -c $(OPTCOMP) fvinject.f


# Pattern rules
%.o : %.f piccom.f fvcom.f;
	$(G77) -c $(OPTCOMP) $*.f

%.o : %.F piccom.f;
	$(G77) -c $(OPTCOMP) $*.F

% : %.f
	$(G77) -o $* $(OPTCOMP) $*.f $(LIB)

% : %.F
	$(G77) -o $* $(OPTCOMP) $*.F $(LIB)


# Distributable archive
sceptic3D.tar.gz : ./accis/libaccisX.a sceptic3D sceptic3Dmpi
	make -C accis mproper
	make -C tools clean
	make clean
	./copyattach.sh
	tar chzf sceptic3D.tar.gz -C .. sceptic3D
	./copyremove.sh


# The following targets will never actually exist
.PHONY: all clean cleandata cleanaccis cleanhdf cleanall ftnchek

all : sceptic3D sceptic3Dmpi

clean :
	-rm *.o
	-rm *.ps
	-rm *.orb
	-rm *.html
	-rm Orbits.txt
	-rm *~

cleandata :
	-rm *.dat
	-rm *.frc
	-rm *.h5

cleanaccis :
	make -C accis clean
	-rm ./accis/libaccisX.a

cleanall :
	make clean
	make cleandata
	make cleanaccis

ftnchek :
	ftnchek -nocheck -nof77 -calltree=text,no-sort -mkhtml -quiet -brief sceptic3D.F *.f

