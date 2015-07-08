#
#  Makefile for isobox 
#
.SUFFIXES: .o .f90 .F90 

include ../Makefile.inc

CUR_DIR = $(CUR_PROJ)/src/pppack
REP_SRC = $(CUR_DIR)

OBJ_MOD = pppack.o 

LIB_OBJECTS = $(OBJ_UTILS) $(OBJ_MOD) 
	  
LIB_OUT = $(REP_LIB)/libpppack.a
	
clean:	
	@( cd $(REP_SRC)/; rm -f *.o *~ *.mod *.a)		

debug:	
	make lib -f Makefile_gfortran F90FLAGS=-g3 FFLAGS=-g3

lib: $(LIB_OBJECTS)
	@ar crv $(LIB_OUT) $(LIB_OBJECTS)
	@cp *.mod $(REP_MOD_TMP)
