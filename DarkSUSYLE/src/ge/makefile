# Makefile for ge directory
# Author: Joakim Edsjo, edsjo@physto.se
# This file is automatically created by makemf.pl on Jun 23, 2006.

# Dependencies and libraries
DINC=../../include

INC_DEP = dsio.h dsidtag.h

vpath %.h $(DINC)

ge = cosd.f dsabsq.f dsbessei0.f dsbessei1.f dsbessek0.f dsbessek1.f \
dsbessek2.f dsbessjw.f dscharadd.f dsf2s.f dsf_int.f dsf_int2.f \
dshiprecint3.f dshunt.f dsi2s.f dsi_trim.f dsidtag.f dsisnan.f \
dsquartic.f dsrnd1.f dsrndlin.f dsrndlog.f dsrndsgn.f dswrite.f erf.f \
erfc.f sind.f spline.f splint.f

all : dsge.o

dsge.o : $(ge) $(INC_DEP) makefile
	cat $(ge) > dsge.f
	$(FF) $(FOPT) -c -I$(DINC) -o dsge.o dsge.f
	rm dsge.f
