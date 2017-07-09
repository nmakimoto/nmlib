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

gtest_io.o: ../include/io.h
gtest_optimization.o: ../include/optimization.h ../include/diff.h
gtest_optimization.o: ../include/matrix.h ../include/random.h
gtest_ode.o: ../include/ode.h ../include/matrix.h
gtest_stat.o: ../include/stat.h ../include/matrix.h
gtest_solver.o: ../include/solver.h ../include/diff.h ../include/matrix.h
gtest_integral.o: ../include/integral.h ../include/random.h
gtest_integral.o: ../include/matrix.h ../include/polynomial.h
gtest_lp.o: ../include/lp.h ../include/matrix.h ../include/io.h
gtest_random.o: ../include/random.h ../include/matrix.h
gtest_matrix_decomp.o: ../include/matrix.h ../include/matrix_decomp.h
gtest_matrix.o: ../include/matrix.h
gtest_sparse.o: ../include/sparse.h ../include/matrix.h
gtest_polynomial.o: ../include/polynomial.h ../include/io.h
gtest_polynomial.o: ../include/solver.h ../include/diff.h ../include/matrix.h
gtest_diff.o: ../include/diff.h ../include/matrix.h
gtest_spline.o: ../include/spline.h
gtest_robot.o: ../include/robot.h ../include/matrix.h
gtest_nmlib.o: ../include/nmlib.h ../include/chrono.h ../include/diff.h
gtest_nmlib.o: ../include/matrix.h ../include/fft.h ../include/integral.h
gtest_nmlib.o: ../include/random.h ../include/polynomial.h ../include/io.h
gtest_nmlib.o: ../include/kalman.h ../include/lp.h ../include/matrix_decomp.h
gtest_nmlib.o: ../include/ode.h ../include/optimization.h ../include/robot.h
gtest_nmlib.o: ../include/route.h ../include/solver.h ../include/sparse.h
gtest_nmlib.o: ../include/spline.h ../include/stat.h
gtest_kalman.o: ../include/kalman.h ../include/diff.h ../include/matrix.h
gtest_fft.o: ../include/fft.h ../include/matrix.h
gtest_route.o: ../include/route.h
