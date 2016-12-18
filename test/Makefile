GTESTDIR=/usr/src/gtest
CC=g++
CXXFLAGS=-Wall -Wextra -O -I../include -I$(GTESTDIR) -std=c++11
LDLIBS=-lpthread
SRC=$(wildcard gtest_*.cpp) random.cpp sparse.cpp gtest-all.cpp
BIN=gtest_main


all : $(BIN)
test: $(BIN)
	./$(BIN)
clean:
	-rm  $(BIN)  gtest-all.cpp  *.o  *.bak  *~
depend:
	makedepend -Y -I. -I../include $(SRC)
$(BIN): $(SRC:.cpp=.o)

gtest-all.cpp:
	ln -s $(GTESTDIR)/src/gtest-all.cc ./gtest-all.cpp


# DO NOT DELETE
