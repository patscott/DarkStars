# Makefile for xcmlib directory
# Author: Joakim Edsjo, edsjo@physto.se
# This file is automatically created by makemf.pl on Jun 23, 2006.

# Dependencies and libraries
DINC=../../include

INC_DEP = 

vpath %.h $(DINC)

xcmlib = d1mach.f dqagse.f dqagseb.f dqelg.f dqk21.f dqk21b.f dqpsrt.f

all : dsxcmlib.o

dsxcmlib.o : $(xcmlib) $(INC_DEP) makefile
	cat $(xcmlib) > dsxcmlib.f
	$(FF) $(FOPT) -c -I$(DINC) -o dsxcmlib.o dsxcmlib.f
	rm dsxcmlib.f
