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


##############################################################################
#  Running collection of things to clean up.
MESS :=


##############################################################################
#  Executables.
BIN     := MD
PROFBIN := MD-prof

MESS += $(BIN) $(PROFBIN)


##############################################################################
#  C source code and generated files.
SRC=MD.c control.c util.c 
OBJ=$(SRC:.c=.o)

MESS += $(OBJ)


##############################################################################
#  C compiler and flags.
CC := pgcc
CFLAGS := -g


##############################################################################
#  Targets.
default: all

all: $(BIN) $(PROFBIN)

$(BIN): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) -lm

$(PROFBIN): $(OBJ)
	$(CC) $(CFLAGS) -pg -o $@ $(OBJ) -lm

output.dat: MD input.dat
	./MD

clean:
	rm -f $(MESS)

$(OBJ) : coord.h Makefile
