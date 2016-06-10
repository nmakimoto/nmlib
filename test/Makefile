GTESTDIR=/usr/src/gtest
CC=g++
CXXFLAGS=-Wall -O -I../include -I$(GTESTDIR)
LDLIBS=-lpthread
SRC=$(wildcard gtest_*.C) random.C sparse.C
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

gtest_diff.o: ../include/diff.H ../include/matrix.H
gtest_fft.o: ../include/fft.H ../include/matrix.H
gtest_io.o: ../include/io.H
gtest_matrix.o: ../include/matrix.H
gtest_matrix_decomp.o: ../include/matrix.H ../include/matrix_decomp.H
gtest_polynomial.o: ../include/polynomial.H ../include/io.H
gtest_polynomial.o: ../include/solver.H ../include/diff.H ../include/matrix.H
gtest_random.o: ../include/random.H ../include/matrix.H ../include/mt64.h
gtest_solver.o: ../include/solver.H ../include/diff.H ../include/matrix.H
gtest_sparse.o: ../include/sparse.H ../include/matrix.H
gtest_spline.o: ../include/spline.H
gtest_stat.o: ../include/stat.H ../include/matrix.H
random.o: mt19937-64.c ../include/mt64.h
sparse.o: ../include/sparse.H ../include/matrix.H
