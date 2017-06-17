GTESTINC=/usr/include
GTESTSRC=/usr/src/gtest

CC=g++
CXXFLAGS=-Wall -Wextra -O -I../include -I$(GTESTINC) -I$(GTESTSRC) -std=c++11
LDLIBS=-lpthread
SRC=$(wildcard gtest_*.cpp) gtest-all.cpp
BIN=gtest_main


all : $(BIN)
test: $(BIN)
	./$(BIN)
clean:
	-rm  $(BIN)  gtest-all.cpp  *.o  *.bak  *~
	#makedepend -Y
depend:
	makedepend -Y -I. -I../include $(SRC)
$(BIN): $(SRC:.cpp=.o)

gtest-all.cpp:
	cp  $(GTESTSRC)/src/gtest-all.cc  ./gtest-all.cpp


# DO NOT DELETE
