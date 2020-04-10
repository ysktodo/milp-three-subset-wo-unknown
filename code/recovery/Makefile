#------------------------------------------------------
# Makefile for C/C++ Program
#------------------------------------------------------
# Target: a.out
#------------------------------------------------------

TARGET=a.out

#------------------------------------------------------
# Default Parameters
#------------------------------------------------------

CC=g++
OPT=-m64 -std=c++17
INC=-I/$$GUROBI_HOME/include/
LIB=-L/$$GUROBI_HOME/lib/ -lgurobi_c++ -lgurobi81 -lm

#------------------------------------------------------
# Compile Option
#------------------------------------------------------

-include makefile.opt

#------------------------------------------------------
# Definition
#------------------------------------------------------

.SUFFIXES:.cpp .c .o .h

#---
# Source Files
#---

SRC=$(shell ls *.cpp)
HED=$(shell ls *.h)
OBJ=$(SRC:.cc=.o)

#------------------------------------------------------
# rules
#------------------------------------------------------

all: $(TARGET)
$(TARGET): $(OBJ)
	$(CC) $(OPT) -o $(TARGET) $(OBJ) $(INC) $(LIB) 

.c.o:
	$(CC) $(OPT) -c $< $(INC)
.cpp.o:
	$(CC) $(OPT) -c $< $(INC)

dep:
	g++ -MM -MG $(SRC) >makefile.depend

clean:
	rm -f $(TARGET) $(TARGET).exe
	rm -f *.o *.obj
	rm -f *~ *.~*

tar:
	tar cvzf $(TARGET).tar.gz $(SRC) $(HED) makefile

#--------------------------------------------------
-include makefile.depend
