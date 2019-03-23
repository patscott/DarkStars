# Makefile for le directory
# Author: Pat Scott, pat@physto.se
# Date: 2007-07-05

# Dependencies and libraries
DINC1=../../include
DINC2=../ini

INC_DEP = dssusy.h dsversion.h dssubversion.h dsdir.h dsidtag.h\
dssun.h dsparam.h

vpath %.h $(DINC1):$(DINC2)

le = dsLEinit.f dsLEfint.f dshiprecint1.f dsntdqagsec.f dsntdqk21c.f \
qromb.f polint.f trapzd.f

nt = dshiprecint.f dshiprecint2.f dsntdqagse.f dsntdqagseb.f dsntdqk21.f dsntdqk21b.f \
dsntsundens.f dsntsunmfrac.f dsntsunread.f dsntsuncdensint.f dsntsuncdfunc.f \
dsntsunpotint.f dsntspfunc.f dsntsunmass.f
nt := $(nt:%=../nt/%)

ge = dshunt.f dsf_int.f spline.f splint.f erf.f
ge := $(ge:%=../ge/%)

xcmlib = d1mach.f dqpsrt.f dqelg.f
xcmlib := $(xcmlib:%=../xcmlib/%)

sources = $(le) $(nt) $(ge) $(xcmlib)

all : dsle.o

dsle.o : $(sources) $(INC_DEP) makefile
	@echo $(INC_DEP) $(nt)
	cat $(sources) > dsle.f
	$(FF) $(FOPT) -c -I$(DINC1) -I$(DINC2) -o dsle.o dsle.f
	#rm dsle.f
