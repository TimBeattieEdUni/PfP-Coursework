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

MESS :=

BIN     := MD
PROFBIN := MD-profiling

MESS += $(BIN) $(PROFBIN)

SRC=MD.c control.c util.c 

OBJ=$(SRC:.c=.o)

MESS += $(OBJ)

CC := pgcc

CFLAGS :=
$(BIN):      CFLAGS := -g
$(PROFBIN):  CFLAGS := -g -pg

default: all

all: $(BIN) $(PROFBIN)

$(BIN): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) -lm

$(PROFBIN): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) -lm


output.dat: MD input.dat
	./MD

clean:
	rm -f $(MESS)

$(OBJ) : coord.h Makefile


