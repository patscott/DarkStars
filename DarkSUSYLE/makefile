# Main makefile for DarkSUSY.
# Technical questions about this makefile should be directed to
# Joakim Edsjo, edsjo@physto.se
# Date: June 23, 2006
# Additions for FeynHiggs by Paolo Gondolo 2002/06/14.

# This file will create libdarksusy.a
# It will also create binary tables from the supplied ascii tables
# to speed up particle yield calculations.
# This file is created from makefile.in by running ./configure.
# Type 'make' to make DarkSUSY
# Type 'make install' to install (default is in /usr/local)

########################################################
### READ THROUGH AND CHANGE THE DEFAULT VALUES BELOW ###
### BEFORE COMPILING.                                ###
########################################################

###### Compiler and compiler options                     ######
###### DarkSUSY root directory and Feynhiggs directories ######

# DS_ROOT is the root directory of DarkSUSY.
# DS_INSTALL is the install directory. When DarkSUSY is installed, you
# can in principle delete the DarkSUSY root directory if you wish.

# If you want the full FeynHiggs in DarkSUSY, set FEYNHIGGSDIR to
# the name of your FeynHiggs directory. This is the directory that 
# contains the source files straight out of the FeynHiggs distribution.
# (You need FeynHiggs version 1.2.2, see http://www.feynhiggs.de)
# If you do not have FeynHiggs, comment out the line setting FEYNHIGGSDIR.

# Note that, on some systems (e.g. Mac OS X), the files in 
# src/xfeynhiggs need the option -femulate-complex, so define
# separate compilation options for these files using the variable FOPTFH. 
# It doesn't seem necessary to compile the full FeynHiggs package with these
# options (the compile is extremeley slow in this case).

# Determine where to install stuff (prefix is set in configure)
prefix=$(PWD)
# DS_INSTALL is where the library and data files will be installed
DS_INSTALL=${prefix}
# DS_ROOT is where the root of the distribution is
DS_ROOT=$(PWD)
#FEYNHIGGSDIR=$(DS_ROOT)/contrib/feynhiggs/
FEYNHIGGSDIR=${prefix}/contrib/feynhiggs

LIB=${DS_INSTALL}/lib

###############################################
#### DO NOT CHANGE ANYTHING BELOW THIS LINE ###
###############################################

darksusy : config_paths darksusy_lib dsfeynhiggs dsgalprop

all : libisajet feynhiggs darksusy galprop

config_paths config :
	perl $(DS_ROOT)/scr/config.pl $(DS_ROOT) $(DS_INSTALL)

darksusy_lib dslib :
	mkdir -p lib
	cd $(DS_ROOT)/src; make all

install : install_builds install_tables

# Only install builds if DS_INSTALL points to some other directory
# than DS_ROOT
install_builds:
	if [ ! -f $(DS_INSTALL)/makefile.in ] ; then \
		echo 'The install directory is different from DS_ROOT' ; \
		echo 'Will install libraries and major data files.' ; \
		mkdir -p $(DS_INSTALL) ; \
		tar cf - lib share include --exclude .svn | (cd $(DS_INSTALL); tar xvf -) ; \
	fi


dsfeynhiggs :
	if [ -d "$(FEYNHIGGSDIR)" ] ; then \
	   (cd $(FEYNHIGGSDIR); make ) ; \
	   cd $(DS_ROOT)/src/xfeynhiggs ; make full ; \
	else \
	   cd $(DS_ROOT)/src/xfeynhiggs ; make fast ; \
	fi

dsgalprop:
	if [ `ar t $(DS_ROOT)/lib/libdarksusy.a dsgalprop.o` ] ; then \
	cd $(DS_ROOT)/src/xgalprop ; make full ; \
	else \
	cd $(DS_ROOT)/src/xgalprop ; make fast ; \
	fi

cfitsio:
	cd $(DS_ROOT)/contrib/cfitsio ; ./configure; make ; \
	cp libcfitsio.a $(DS_INSTALL)/lib/ ; \

really_all : all libisajet

install_tables :
	mkdir -p $(DS_INSTALL)/share/DarkSUSY
	cd $(DS_ROOT)/share/DarkSUSY; make all

libisajet :
	cd $(DS_ROOT)/contrib/isajet; make

feynhiggs :
	cd $(DS_ROOT)/contrib/feynhiggs; make

galprop:
	cd $(DS_ROOT)/src/xgalprop ; make full
#	cd $(DS_ROOT)/contrib/galprop/v42.3ds; make
	cd $(DS_ROOT)/contrib/galprop/v50p; make
	cd $(DS_ROOT); if [ -f .no_cfitsio ]; then\
		make cfitsio; rm .no_cfitsio; fi

clean :
	cd $(DS_ROOT)/src; $(MAKE) clean
	cd $(DS_ROOT)/share/DarkSUSY; $(MAKE) clean
	cd $(DS_ROOT)/docs; rm -f *.aux; rm -f *.log; rm -f *.dvi
	cd $(DS_ROOT)/contrib/isajet; $(MAKE) clean
#	cd $(DS_ROOT)/contrib/galprop/v42.3ds; $(MAKE) clean
	cd $(DS_ROOT)/contrib/galprop/v50p; $(MAKE) clean
	cd $(DS_ROOT)/contrib/feynhiggs; $(MAKE) new

clobber distclean : clean rmmore

# rmmore removes things added in the DarkSUSY root during configure
# and build
rmmore:
	rm -f $(DS_ROOT)/lib/libdarksusy.a
	rm -f $(DS_ROOT)/lib/aldata.o
	rm -f $(DS_ROOT)/lib/libcfitsio.a
	rm -f $(DS_ROOT)/share/DarkSUSY/ascii2bin
	rm -f $(DS_ROOT)/share/DarkSUSY/sim*.bin
	find . -name '#*' -exec rm -f {} \;
	find . -name '*~' -exec rm -f {} \;
	find . -name '._*' -exec rm -f {} \;
	find . -name '*.o' -exec rm -f {} \;
	find test -perm +u=x -type f -exec rm -f {} \;
	rm -f `find . -name makefile.in | sed s/\.in//`
	rm -f config.log config.status
	rm -fr autom4te.cache
	rm -f src/ini/dsversion.h
	rm -f src/ini/dsdir.h
	rm -f contrib/galprop/v42.3ds/Configure.cc
	rm -f contrib/galprop/v42.3ds/nuc_package.cc
	rm -f contrib/galprop/v42.3ds/WNEWTR_FUNC_aws.f
	rm -f .no_cfitsio
	cd contrib/cfitsio; make -k distclean; cd $(DS_ROOT)
#	find . -name '.[a-z,A-Z,#]*' -exec rm -f {} \;

# uninstall removes things in the install directory
# Note, we have to check if the install directory is the current working
# directory, in which case some things should be kept
# We do this by checking if makefile.in is in the install directory.
# If it is, we know that DS_INSTALL must point to DS_ROOT.
uninstall:
	rm -f $(DS_INSTALL)/lib/libdarksusy.a
	rm -f $(DS_INSTALL)/lib/aldata.o
	if [ ! -f $(DS_INSTALL)/makefile.in ] ; then \
		echo 'The install directory is different from DS_ROOT' ; \
		echo 'Will clean up the install directory.' ; \
		rm -f `find include -name '*.h' | sed s#^#${DS_INSTALL}/#` ; \
		rm -fr $(DS_INSTALL)/share/DarkSUSY ; \
	fi

docs: pdf-manual html-manual

pdf-manual pdf-Manual :
	scr/headers2tex.pl
	cd $(DS_ROOT)/docs; latex Manual; latex Manual; latex Manual
	cd $(DS_ROOT)/docs; dvips -o Manual.ps Manual
	cd $(DS_ROOT)/docs; ps2pdf Manual.ps

pdf-manual-short pdf-Manual-short :
	scr/headers2tex.pl --no-headers
	cd $(DS_ROOT)/docs; \
	latex Manual-short; latex Manual-short; latex Manual-short
	cd $(DS_ROOT)/docs; dvips -o Manual-short.ps Manual-short
	cd $(DS_ROOT)/docs; ps2pdf Manual-short.ps

html-manual :
	scr/headers2tex.pl
	cd $(DS_ROOT)/docs; latex2html -white Manual.tex

.PHONY: all config_paths darksusy_lib dslib private install_tables \
	clean clobber uninstall rmmore libisajet pdf-manual html-manual

# Inserted for DarkStars; Pat Scott 28-06-2007
darksusyLE_lib :
	mkdir -p lib
	cd $(DS_ROOT)/src; make libdarksusyLE.a

cleanLE : 
	cd $(DS_ROOT)/src; make cleanLE
# End DarkStars section
