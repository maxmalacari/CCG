##################################################
# makefile for the Coverage & Anisotropy Toolkit #
##################################################

SYSTEM=$(shell uname)
CC = gcc 
CXX = g++
OPTFLAGS = #-O3 -march=native
DBGFLAGS = -g3
WFLAGS = -D__USE_FIXED_PROTOTYPES__ -Wall
OBJ = ./

.PHONY: clean
#------------------- defs ------------------------------
SRC_DIR = /Users/max/Documents/Work/Anisotropy/Starburst/Toolkit-collaboration-2017-03-15/src/

CFITSIO_INCDIR = $(CFITSIO_DIR)/include
CFITSIO_LIBDIR = $(CFITSIO_DIR)/lib

CHEALPIX_INCDIR = $(HEALPIX_DIR)/include
CHEALPIX_LIBDIR = $(HEALPIX_DIR)/lib

INCDIR = -I$(CFITSIO_INCDIR) \
	-I$(CHEALPIX_INCDIR) \
	-I$(SRC_DIR) \
	$(shell root-config --cflags) \
	-I$(HEALPIX_DIR)/src/cxx/osx/include

LIBDIR = -L${CHEALPIX_LIBDIR} -lchealpix -lm \
	-L${CFITSIO_LIBDIR} -lcfitsio -lMinuit \
	$(shell root-config --glibs) \
	-L$(HEALPIX_DIR)/src/cxx/osx/lib -lhealpix_cxx -lcxxsupport
#-------------------------------------------------------

#------- alias -----------------------------------------
execs = \
	simulateCCG
exeobjs = $(patsubst %.exe,%.o,$(execs))

HEADERS = $(patsubst %.o,%.h,$(libobjs))

thelib = $(SRC_DIR)/libtkit.a
#-------------------------------------------------------

#-------- rules ----------------------------------------
# rules for the executable sources
$(OBJ)%.o:%.cc
	$(COMPILE.cc) $(DBGFLAGS) $(OPTFLAGS) $(WFLAGS) $(INCDIR) -o $@ $<
#-------------------------------------------------------

#------- targets ---------------------------------------
all :	lib $(execs) 
lib :
	cd $(SRC_DIR)/; make
clean :
	cd $(SRC_DIR)/; make clean;
	@echo "Deleting library objects, executables and associated objects."
	@/bin/rm -f $(thelib) $(execs) $(libobjs) $(exeobjs)
	@/bin/rm -f *~
	@echo "Done"
#-------------------------------------------------------

#-------- specific rules -------------------------------
simulateCCG : simulateCCG.o $(thelib)
	$(CXX) -o $@ $^ $(LIBDIR)
#-------------------------------------------------------

