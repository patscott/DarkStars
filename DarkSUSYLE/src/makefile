# Makefile for the DarkSUSYLE library
#
# Pat Scott, pat@fysik.su.se

LIB=../lib
INC=../include

vpath %.a $(LIB)

DSLEDIRS = cap le

DSLEOBJ = $(DSLEDIRS:%=%_obj)

DSLEOBJWDIR = cap/dscap.o le/dsle.o

$(DSLEOBJ) :
	cd $(@:_obj=) ; make all

libdarksusyLE.a : $(DSLEOBJ)
	ar rv $(LIB)/libdarksusyLE.a $(DSLEOBJWDIR)
	ranlib $(LIB)/libdarksusyLE.a
	
cleanlibLE :
	-rm $(LIB)/libdarksusyLE.a

cleanLE : clean cleanlibLE
	
clean : 
	rm -f */*.o 
