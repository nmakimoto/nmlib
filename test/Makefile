GTESTDIR=/usr/src/gtest
CC=g++
CXXFLAGS=-Wall -Wextra -O -I../include -I$(GTESTDIR) -std=c++11
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

gtest_io.o: ../include/io.H
gtest_lp.o: ../include/lp.H ../include/matrix.H ../include/io.H
gtest_nmlib.o: ../include/nmlib.H ../include/diff.H ../include/matrix.H
gtest_nmlib.o: ../include/fft.H ../include/io.H ../include/lp.H
gtest_nmlib.o: ../include/matrix_decomp.H ../include/polynomial.H
gtest_nmlib.o: ../include/random.H ../include/robot.H ../include/solver.H
gtest_nmlib.o: ../include/sparse.H ../include/spline.H ../include/stat.H
gtest_fft.o: ../include/fft.H ../include/matrix.H
gtest_solver.o: ../include/solver.H ../include/diff.H ../include/matrix.H
gtest_random.o: ../include/random.H ../include/matrix.H
gtest_sparse.o: ../include/sparse.H ../include/matrix.H
gtest_matrix.o: ../include/matrix.H
gtest_polynomial.o: ../include/polynomial.H ../include/io.H
gtest_polynomial.o: ../include/solver.H ../include/diff.H ../include/matrix.H
gtest_diff.o: ../include/diff.H ../include/matrix.H
gtest_spline.o: ../include/spline.H
gtest_robot.o: ../include/robot.H ../include/matrix.H ../include/io.H
gtest_stat.o: ../include/stat.H ../include/matrix.H
gtest_matrix_decomp.o: ../include/matrix.H ../include/matrix_decomp.H
random.o: ../include/random.H ../include/matrix.H
sparse.o: ../include/sparse.H ../include/matrix.H
