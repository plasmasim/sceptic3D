# Universal Makefile for sceptic3D


# Set shell to bash (default is sh)
SHELL := /bin/bash

# Set number of threads to use for HDF make
NUMPROC := 8

# Set compilers
G77 := $(shell ./setcomp f77)
G77nonmpi := $(shell ./setcomp f77 nonmpi)
G90 := $(shell ./setcomp f90)
G90nonmpi := $(shell ./setcomp f90 nonmpi)

# HDF root directory
#HDFDIR := $(realpath hdf5-1.8.4)
# realpath not available on loki, so use hack
DIRHDF := $(PWD)/hdf5-1.8.4

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

# Libraries and options to pass to linker for HDF version
LIBHDF := $(LIB)
# HDF libraries and options (from h5fc -show in DIRHDF/bin)
LIBHDF += -L$(DIRHDF)/lib -lhdf5hl_fortran -lhdf5_hl \
          -lhdf5_fortran -lhdf5 -lz -lm -Wl,-rpath -Wl,$(DIRHDF)/lib
# Note that the -Wl,-rpath... options are needed since LD_LIBRARY_PATH
#   does not contain the location of the hdf shared libraries at runtime 


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

# Options to pass to compiler for HDF version
OPTCOMPHDF := $(OPTCOMP)
# Include directory with HDF modules
OPTCOMPHDF += -I$(DIRHDF)/include
# Enable HDF by defining 'HDF' for pre-compiler
OPTCOMPHDF += -DHDF

# Options to pass to compiler for MPI version
OPTCOMPMPI := $(OPTCOMP)
# Enable MPI by defining 'MPI' for pre-compiler
OPTCOMPMPI += -DMPI

# Options to pass to compiler for MPI & HDF version
OPTCOMPMPIHDF := $(OPTCOMPHDF)
# Enable MPI by defining 'MPI' for pre-compiler
OPTCOMPMPIHDF += -DMPI


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

# Objects for HDF version of sceptic3D
OBJHDF := $(OBJ) \
          outputhdf.o

# Objects for MPI version of sceptic3D
OBJMPI := $(OBJ) \
          cg3dmpi.o \
          mpibbdy.o \
          shielding3D_par.o

# Objects for MPI & HDF version of sceptic3D
OBJMPIHDF := $(OBJMPI) \
          outputhdf.o

# Default target is serial sceptic3D without HDF support
sceptic3D : sceptic3D.F piccom.f $(OBJ) ./accis/libaccisX.a
	$(G77) $(OPTCOMP) -o sceptic3D sceptic3D.F $(OBJ) $(LIB)

# sceptic3D with HDF
sceptic3Dhdf : sceptic3D.F piccom.f $(OBJHDF) ./accis/libaccisX.a
	$(G77) $(OPTCOMPHDF) -o sceptic3Dhdf sceptic3D.F $(OBJHDF) $(LIBHDF)

# sceptic3D with MPI
sceptic3Dmpi : sceptic3D.F piccom.f piccomcg.f $(OBJMPI) ./accis/libaccisX.a
	$(G77) $(OPTCOMPMPI) -o sceptic3Dmpi sceptic3D.F $(OBJMPI) $(LIB)

# sceptic3D with MPI & HDF
sceptic3Dmpihdf : sceptic3D.F piccom.f piccomcg.f $(OBJMPIHDF) ./accis/libaccisX.a
	$(G77) $(OPTCOMPMPIHDF) -o sceptic3Dmpihdf sceptic3D.F $(OBJMPIHDF) $(LIBHDF)


# HDF related rules
outputhdf.o : outputhdf.f piccom.f colncom.f $(DIRHDF)/lib/libhdf5.a
	$(G90) -c $(OPTCOMPHDF) outputhdf.f

# Though more than one hdf library used, choose one as trigger
$(DIRHDF)/lib/libhdf5.a :
	cd $(DIRHDF) &&	\
	./configure --prefix=$(DIRHDF) --enable-fortran \
	FC=$(G90nonmpi) && \
	make -j$(NUMPROC) && \
	make install
# Note that providing an mpi compiler to hdf will cause it to build
#   the MPI version, which is not needed (and didn't work on sceptic)


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

all : sceptic3D sceptic3Dhdf sceptic3Dmpi sceptic3Dmpihdf

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

cleanhdf :
	make -C $(DIRHDF) clean
	-rm $(DIRHDF)/lib/libhdf5.a

cleanall :
	make clean
	make cleandata
	make cleanaccis
	make cleanhdf

ftnchek :
	ftnchek -nocheck -nof77 -calltree=text,no-sort -mkhtml -quiet -brief sceptic3D.F *.f

