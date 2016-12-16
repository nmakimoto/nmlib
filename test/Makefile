GTESTDIR=/usr/src/gtest
CC=g++
CXXFLAGS=-Wall -Wextra -O -I../include -I$(GTESTDIR) -std=c++11
LDLIBS=-lpthread
SRC=$(wildcard gtest_*.cpp) random.cpp sparse.cpp
BIN=gtest_main


all : $(BIN)
test: $(BIN)
	./$(BIN)
clean:
	-rm  $(BIN)  *.o  *.bak  *~
depend:
	makedepend -Y -I. -I../include $(SRC)
$(BIN): $(SRC:.cpp=.o)


# DO NOT DELETE
