GTESTDIR=/usr/src/gtest
CC=g++
CXXFLAGS=-Wall -O -I../include -I$(GTESTDIR)
SRC=$(wildcard gtest_*.C)
BIN=gtest_main


all : $(BIN)
test: $(BIN)
	./$(BIN)
clean:
	-rm  $(BIN)  *.o  *.bak  *~
depend:
	makedepend -Y -I. -I../include $(SRC)
$(BIN): $(SRC:.C=.o)


# DO NOT DELETE
