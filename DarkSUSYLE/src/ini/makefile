# Makefile for ini directory
# Author: Joakim Edsjo, edsjo@physto.se
# This file is automatically created by makemf.pl on Jun 23, 2006.

# Dependencies and libraries
DINC=../../include

INC_DEP = dsprep.h dsrncom.h dsascom.h dsaccom.h dsandwcom.h dsdir.h \
dsidtag.h dspbcom.h dsversion.h dssusy.h

vpath %.h $(DINC)

ini = dscval.f dsfval.f dsinit.f dsival.f dskillsp.f dslowcase.f \
dslval.f dsreadpar.f

all : dsini.o

dsini.o : $(ini) $(INC_DEP) makefile
	cat $(ini) > dsini.f
	$(FF) $(FOPT) -c -I$(DINC) -o dsini.o dsini.f
	rm dsini.f
