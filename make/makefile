# DarkStars makefile
#
# This is a heavily modified version of the original EZ makefile by Bill Paxton.
#
# Pat Scott, pat@fysik.su.se

# NOTE on fixed versus free format.
# I've written the code so it works either way (I hope).
# For free format, you don't need to do anything about line length.
# But if you want to use fixed format, it does need 132 column lines.
# Where possible, I've already put in compiler flags to allow long lines.
# If you have no idea what any of this is about, just try compiling.
# Then if you get some syntax errors in long lines, you'll need to find the right flags.

# If your compiler is listed below, just delete the "#" in the FC line for it.

# for the Absoft compiler
#FC = f95 -W132

# for the Intel compiler
#FC = ifort -w -132
# -w disables warnings
# -132 allows long lines

# for the Intel compiler with flags set for Intel Core Duo processor (laptop)
# FC = ifort -132 -xP -arch SSE2 -vec-report0
# for the Intel compiler with flags set for Intel Pentium 4 with no SSE extensions (GLAST boxes)
# FC = ifort -132 -vec-report0 -xW
# for the Intel compiler with flags set for Intel Core 2 Duo / Xeon processors (OSX 8 core machine)
# FC = ifort -132 -vec-report0 -xT
# -axP generates code for Intel core duo processors.
# -xP generates code for Intel core duo processors.
# -ip enables additional interprocedural optimizations
# -p for grof profiling
# -ipo for extra interfile optimization and inlining
# -vec-report0 to eliminate diagnostics about vectorized loops

# for gfortran
FC = gfortran -ffixed-line-length-132

# for the IBM compiler with free format
#FC = xlf95
# for the IBM compiler with fixed format
#FC = xlf95 -qfixed=132

# for the Sun compiler.  best of luck.  so far no one has been able to make it work.

# for the Lahey compiler
#FC = lf95 --wide (as of Aug 14, 2004, this compiles and loads, but seems to give unreliable results.)

# for the g95 compiler (as of Sept 3, 2004, this is working -- but the resulting EZ is slower by more than a factor of 2)
#FC = g95 -ffixed-line-length-132

# use standard optimization
FFLAGS = -O2 #-check all -traceback -fpstkchk -std90 -warn unused -warn declarations
#FFLAGS = -O2 -Rb -Rs
# -Rb checks array boundaries, -Rs checks array conformance
# for performance profiling and gprof with Absoft compiler (others use -pg)
#FFLAGS = -P -g
# then do gprof EZ > profile.doc; open profile.doc

FF = $(FC)
export FC FF FFLAGS

PROG = EZ
RUNDIR = ../run
SRCDIR = ../src
LADIR = ../LAPACK

OBJ_EARLY = star_constants.o star_controls.o star_data.o ez_cycle_data.o ez_magnitude_data.o star_extras.o
OBJ_LATE = \
 ez_vcool_data.o ez_shell_data.o ez_opacity_data.o ez_ionization_data.o \
 ez_state_data.o ez_nuclear_data.o ez_equations_data.o ez_solve_data.o ez_do_one_data.o ez_data.o \
 ez_utils.o ez_log_utils.o ez_log.o \
 ez_opacity.o ez_convection.o ez_nuclear.o ez_vcool.o ez_ionization.o ez_state.o \
 ez_shell.o ez_equations.o ez_magnitude.o ez_report.o ez_solve.o ez_cycle.o ez_setup.o \
 ez_flash.o ez_driver.o ez_do_one_utils.o
OBJ_LATER = ez_do_one.o ez_example.o interp_1D.o demo1.o demo2.o experiment.o
OBJ_VERY_LATE = dispatch.o
OBJ = $(OBJ_EARLY) $(OBJ_LATE) $(OBJ_LATER) $(OBJ_VERY_LATE)

# the following LAPACK files are used for doing the array solution

LAOBJ = \
 dcopy.o lsame.o dgemm.o dgemv.o dger.o dscal.o dswap.o dtbsv.o dtrsm.o idamax.o xerbla.o \
 ieeeck.o ilaenv.o dgbsv.o dgbtf2.o dgbtrf.o dgbtrs.o dlaswp.o \
 dasum.o daxpy.o ddot.o dgbcon.o dgbequ.o dgbmv.o dgbrfs.o dgbsvx.o dlassq.o \
 dlabad.o dlacn2.o dlacpy.o dlamch.o dlangb.o dlantb.o dlaqgb.o dlatbs.o drscl.o iparmq.o

# Inserted for DarkStars; Pat Scott July 2007

DS_ROOT = ../DarkSUSYLE
DS_INSTALL = $(DS_ROOT)
DS_INC = $(DS_ROOT)/include
DS_LIB = $(DS_INSTALL)/lib
DSEDITION = darksusyLE

NUMDIR = ../num
# The commented routines are derived from proprietry code and therefore not included in the public DarkStars release.
# You can replace them with your own personal versions, or just rely on Simpson's rule for all integrals
# (which is the default anyway).  RK45 is faster, but not as reliable as Simpson's rule for this problem.
NUMOBJ = zeroin.o #zbrent.o #fehl.o rkfs.o rkf45.o DkStrs_Romberg.o

TSDIR = ../TSPACK
TSOBJ = ENDSLP.o SIGS.o SNHCSH.o STORE.o YPCOEF.o YPC1.o YPC1P.o YPC2.o YPC2P.o TSPSI.o \
 INTRVL.o HVAL.o HPVAL.o

DSTARSPROG = DarkStars
DSTARSOBJ_VERY_EARLY = DkStrs_akraw.o DkStrs_constants.o DkStrs_data.o
DSTARSOBJ_EARLY = DkStrs_utils.o \
 DkStrs_WIMPdens.o DkStrs_annihilation.o DkStrs_fint_ext.o DkStrs_transport.o
DSTARSOBJ_MIDDLE = DkStrs_init.o DkStrs_admin.o
DSTARSOBJ_LATE  = DarkStars.o
DSTARSOBJ = $(DSTARSOBJ_VERY_EARLY) $(DSTARSOBJ_EARLY) $(DSTARSOBJ_MIDDLE) $(DSTARSOBJ_LATE)

all :  darkstars

darkstars : $(LAOBJ) $(TSOBJ) $(DSEDITION:%=%_lib) $(OBJ_EARLY) $(DSTARSOBJ_VERY_EARLY) $(NUMOBJ) \
	$(DSTARSOBJ_EARLY) $(OBJ_LATE) $(DSTARSOBJ_MIDDLE) $(OBJ_LATER) $(DSTARSOBJ_LATE) $(OBJ_VERY_LATE)
	$(FC) $(FFLAGS) -L$(DS_LIB) -o $(RUNDIR)/$(DSTARSPROG) \
	$(LAOBJ) $(TSOBJ) $(NUMOBJ) $(OBJ) $(DSTARSOBJ) -l$(DSEDITION)

$(DSEDITION:%=%_lib):
	cd $(DS_ROOT); make $(DSEDITION:%=%_lib)

$(DSTARSOBJ) : %.o : $(SRCDIR)/%.f $(DSEDITION:%=%_lib)
	$(FC) -c $(FFLAGS) -I$(DS_INC) -L$(DS_LIB) -l$(DSEDITION) $<

$(NUMOBJ) : %.o : $(NUMDIR)/%.f
	$(FC) $(FFLAGS) -c $<

$(TSOBJ) : %.o : $(TSDIR)/%.f
	$(FC) $(FFLAGS) -c $<

clean : clean_EZ
	-@rm -f $(DSTARSOBJ) $(TSOBJ) $(NUMOBJ)

cleanall : cleanall_EZ clean_darksusy clean

clean_darksusy :
	cd $(DS_ROOT); make cleanLE

# End DarkStars section

%.o: $(SRCDIR)/%.f
	$(FC) -c $(FFLAGS) $<

%.o: $(SRCDIR)/%.f90
	$(FC) -c $(FFLAGS) $<

%.o: $(LADIR)/%.f
	$(FC) -c $(FFLAGS) $<

clean_EZ :
	-@rm -f $(OBJ) *.mod

cleanall_EZ :
	-@rm -f $(LAOBJ) $(OBJ) *.mod
