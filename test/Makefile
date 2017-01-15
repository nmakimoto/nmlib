GTESTINC=/usr/include
GTESTSRC=/usr/src/gtest

CC=g++
CXXFLAGS=-Wall -Wextra -O -I../include -I$(GTESTINC) -I$(GTESTSRC) -std=c++11
LDLIBS=-lpthread
SRC=$(wildcard gtest_*.cpp) random.cpp sparse.cpp gtest-all.cpp
BIN=gtest_main


all : $(BIN)
test: $(BIN)
	./$(BIN)
clean:
	-rm  $(BIN)  gtest-all.cpp  random.cpp  sparse.cpp  *.o  *.bak  *~
	#makedepend -Y
depend:
	makedepend -Y -I. -I../include $(SRC)
$(BIN): $(SRC:.cpp=.o)

gtest-all.cpp:
	cp  $(GTESTSRC)/src/gtest-all.cc  ./gtest-all.cpp
random.cpp: ../src/random.cpp
	cp  ../src/random.cpp ./
sparse.cpp: ../src/sparse.cpp
	cp  ../src/sparse.cpp ./


# DO NOT DELETE
