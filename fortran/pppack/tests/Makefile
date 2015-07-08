#
#  Makefile for geometries/tests
#
.SUFFIXES: .o .f90 .F90 

include ../../Makefile.inc

CUR_DIR = $(CUR_PROJ)/src/pppack/tests
REP_SRC = $(CUR_DIR)
REP_BIN = $(CUR_DIR)

OBJECTS = tests.o main.o

INCLUDE := $(INCLUDE) -I$(REP_LIB) -I$(REP_MOD) -I$(REP_MOD_TMP)
LIB_DIR := $(LIB_DIR) -L$(REP_LIB)
LIB := $(LIB)   -lpppack -ltracelog
	
all:	$(OBJECTS)
	$(LD) $(F90FLAGS) $(OBJECTS) $(LIB_DIR)  $(LIB) -o $(REP_BIN)/run

clean:	
	@( cd $(REP_SRC)/; rm -f *.o *~ *.mod *.a)		

debug:	
	make lib -f Makefile_gfortran F90FLAGS=-g3 FFLAGS=-g3

