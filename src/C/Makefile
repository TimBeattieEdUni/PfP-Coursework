##############################################################################
# MSC HPC
# Performance Programming
# Coursework
# Exam no. B063622
#
# Makefile for MD program.
#
# Targets:
#   default:  MD
#   MD        Builds the MD program.
#


##############################################################################
# Running collection of things to clean up.
#
MESS :=


##############################################################################
# Executables.
#
BIN     := MD

MESS += $(BIN)


##############################################################################
# C source code and generated files.
#
SRC=MD.c control.c util.c
OBJ=$(SRC:.c=.o)

MESS += $(OBJ)


##############################################################################
# C compiler and flags.
#
CC := pgcc
CFLAGS := -fast -Mipa=fast,inline

$(PROFBIN): CFLAGS += -pg


##############################################################################
# Targets.  Note: "make clean" is required when switching between profiling
# and regular builds.
#
default: $(BIN)

$(BIN): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) -lm

output.dat: MD input.dat
	./MD

clean:
	rm -f $(MESS)

$(OBJ) : coord.h Makefile
