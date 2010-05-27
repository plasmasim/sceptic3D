# Universal Makefile for sceptic3D

# Set shell to bash (default is sh)
SHELL := /bin/bash

# Set compilers
G77 := $(shell ./setcomp f77)
G90 := $(shell ./setcomp f90)
G90nonmpi := $(shell ./setcomp f90 nonmpi)

# Set Xlib location
XLIB := $(shell ./setxlib)

# Accis lib location
ACCISLIB ?= ./accis

LIBRARIES :=  -L$(XLIB) -L$(ACCISLIB) -laccisX -lXt -lX11 

# Location of hdf5
HDFDIR := $(realpath hdf5-1.8.4)
# To figure out what to use for the hdf includes and libraries
# run the h5fc script with -show ($(HDFDIR)/bin/h5fc)
HDFINCLUDE := -I$(HDFDIR)/include
HDFLIBRARIES := -L$(HDFDIR)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5

#Default No Warnings
ifeq ("$(NOWARN)","")
	NOWARN=
endif

COMPILE-SWITCHES =-Wall -Wno-unused-variable  $(NOWARN)  -O2  -I.
## For debugging, don't use optimization
#COMPILE-SWITCHES =-Wall -Wno-unused-variable  $(NOWARN)  -I.
# For debugging.
#  -g  -ffortran-bounds-check
# For profiling
#COMPILE-SWITCHES = -Wall -O2 -pg

REINJECT=orbitinject.o extint.o maxreinject.o mcreinject.o

MPICOMPILE-SWITCHES = -DMPI $(COMPILE-SWITCHES)

OBJECTS = initiate.o advancing.o randc.o randf.o diags.o outputs.o outputhdf.o	\
 chargefield.o $(REINJECT) stringsnames.o		\
rhoinfcalc.o shielding3D.o utils.o


MPIOBJECTS=cg3dmpi.o mpibbdy.o shielding3D_par.o

all : sceptic3D

sceptic3D :  sceptic3D.F  piccom.f errcom.f  ./accis/libaccisX.a $(OBJECTS)
	$(G77) $(COMPILE-SWITCHES) $(HDFINCLUDE) $(HDFLIBRARIES) -o sceptic3D sceptic3D.F  $(OBJECTS) $(LIBRARIES)

sceptic3Dmpi : sceptic3D.F  piccom.f errcom.f piccomcg.f ./accis/libaccisX.a $(OBJECTS) $(MPIOBJECTS)
	$(G77) $(MPICOMPILE-SWITCHES) $(HDFINCLUDE) $(HDFLIBRARIES) -o sceptic3Dmpi  sceptic3D.F   $(OBJECTS) $(MPIOBJECTS) $(LIBRARIES)

./accis/libaccisX.a : ./accis/*.f
	make -C accis

orbitint : orbitint.f coulflux.o $(OBJECTS) ./accis/libaccisX.a
	$(G77) $(COMPILE-SWITCHES) -o orbitint orbitint.f $(OBJECTS) coulflux.o $(LIBRARIES)

coulflux.o : tools/coulflux.f
	$(G77) -c $(COMPILE-SWITCHES) tools/coulflux.f

fvinjecttest : fvinjecttest.F fvinject.o reinject.o initiate.o advancing.o chargefield.o randf.o fvcom.f
	$(G77)  -o fvinjecttest $(COMPILE-SWITCHES) fvinjecttest.F fvinject.o reinject.o initiate.o advancing.o chargefield.o randf.o  $(LIBRARIES)

fvinject.o : fvinject.f fvcom.f piccom.f errcom.f
	$(G77) -c $(COMPILE-SWITCHES) fvinject.f

outputhdf.o : outputhdf.f piccom.f errcom.f colncom.f hdf
	$(G90) -c $(COMPILE-SWITCHES)  $(HDFINCLUDE) outputhdf.f

#pattern rule
%.o : %.f piccom.f errcom.f fvcom.f;
	$(G77) -c $(COMPILE-SWITCHES) $*.f

%.o : %.F piccom.f errcom.f;
	$(G77) -c $(COMPILE-SWITCHES) $*.F

% : %.f
	$(G77)  -o $* $(COMPILE-SWITCHES) $*.f  $(LIBRARIES)

% : %.F
	$(G77)  -o $* $(COMPILE-SWITCHES) $*.F  $(LIBRARIES)

sceptic3D.tar.gz : ./accis/libaccisX.a sceptic3D sceptic3Dmpi
	make -C accis mproper
	make -C tools clean
	make clean
	./copyattach.sh
	tar chzf sceptic3D.tar.gz -C .. sceptic3D
	./copyremove.sh


.PHONY: all clean ftnchek

clean :
	rm -f *.o
	rm -f *.ps
	rm -f *.orb
	rm -f *.html
	rm -f Orbits.txt
	rm -f *~

cleanall :
	make clean
	rm -f *.dat
	rm -f *.frc

ftnchek :
	ftnchek -nocheck -nof77 -calltree=text,no-sort -mkhtml -quiet -brief sceptic3D.F *.f

$(HDFDIR)/lib/libhdf5.a :
	cd $(HDFDIR) &&	\
	./configure --prefix=$(PWD) --enable-fortran \
	FC=$(G90nonmpi) && \
	make -j8 && \
	make install
