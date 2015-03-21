##############################################################################
# MSC HPC
# Performance Programming
# Coursework
# Exam no. B063622
#
# Makefile for MD program.
#
# Targets:
#    default:  MD
#    MD        Builds the MD program.
#    perf      Builds several versions with varying compiler flags and creates 
#              graphs of results.  (not implemented yet)
#

SRC=MD.c control.c util.c 
OBJ=$(SRC:.c=.o)
CC=pgcc


CFLAGS= -g -pg


all: MD

MD: $(OBJ)
	$(CC) $(CFLAGS)  -o $@  $(OBJ) -lm


output.dat: MD input.dat
	./MD


clean:
	rm -f MD $(OBJ) 

$(OBJ) : coord.h Makefile


