#'make' build executable file 'corr'
# 'make clean' remove all .o and executable files
#

IDIR=../include
CC=g++
CFLAGS=-I$(IDIR) -Wall -g
#-Wall is for "warn all". Turns on all warnings g++ can tell you about.
#-g turns on information mode for gdb compiler.

ODIR=obj
#LDIR = ../lib

LIBS=-lhealpix_cxx

_DEPS = corrfunc.h helper.h
DEPS = $(patsubst %, $(IDIR)/%, $(_DEPS))

_OBJ = main.o corrfunc.o helper.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))


#define the executable file
EDIR = ../bin
_MAIN = corr
MAIN = $(patsubst %, $(EDIR)/%, $(_MAIN))

#generic makefile
#--------------------

.PHONY: clean

all : $(OBJ)
	$(CC) -o $(MAIN) $^ $(CFLAGS) $(LIBS)

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(FLAGS)
clean:
	rm -f $(ODIR)/*.o $(MAIN)






