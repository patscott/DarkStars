# Makefile for cap directory
# Author: Pat Scott, pat@physto.se

# Dependencies and libraries
DINC=../../include

INC_DEP = dsparam.h dscapstar.h

vpath %.h $(DINC)

cap = dscapcsint.f dscapcsint1.f dscapcsint2.f dscapfoveru.f \
dscapsetup.f dscapspfunc.f dscapstar.f dscapstardenscomp.f \
dscapstarpotcalc.f dscapstari.f dscapstarmass.f dscapstarpot.f \
dscapstarpotint.f dscapstarvesc.f dscapudfgauss_star.f \
dscapMBhydrogen.f dscapMBheavy.f

all : dscap.o

dscap.o : $(cap) $(INC_DEP) makefile
	cat $(cap) > dscap.f
	$(FF) $(FOPT) -c -I$(DINC) -o dscap.o dscap.f
	rm dscap.f
