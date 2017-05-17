#CC = g++
#CC = xlcpp_r
CC = icpc

#MIC = icpc -mmic

CFLAGS = -xHOST -O2

LDFLAGS = -qopenmp
#LDFLAGS = -fopenmp


#MP = g++
#MP = mpiicpc
MP = mpiicpc
#MP = mpixlcxx
#MP = /usr/lib64/openmpi/bin/mpicxx

headers =	epf.h matrix.h data.h opt.h solve.h
CSOURCES = objec.cpp tauchen.cpp matrixu.cpp epf.cpp regress.cpp utilities.cpp inference.cpp

intermediateC = $(CSOURCES:.cpp=.o)

.cpp.o: $(headers)
	$(CC) ${CFLAGS} $(LDFLAGS) $< -c

matrix-fast :
		cp matrix_fast.h matrix.h
		cp matrixu_fast.cpp matrixu.cpp

matrix-slow :
		cp matrix_slow.h matrix.h
		cp matrixu_slow.cpp matrixu.cpp

monte.o : ${monte} $(headers)
	${MP} ${CFLAGS} ${LDFLAGS} -c ${monte} ${LDFLAGS} -o monte.o

estim.o : ${estim} $(headers)
	${CC} ${CFLAGS} ${LDFLAGS} -c ${estim} ${LDFLAGS} -o estim.o

trial.o : trial.cpp $(headers)
	${CC} ${CFLAGS} ${LDFLAGS} -c trial.cpp ${LDFLAGS}

optim.o : ${optim} opt.h $(headers)
	${CC} ${CFLAGS} ${LDFLAGS} -c ${optim} ${LDFLAGS} -o optim.o

model.o : $(model) $(headers)
	${CC} ${CFLAGS} ${LDFLAGS} -c $(model) ${LDFLAGS} -o model.o

simulate.o : ${simulate} ${headers}
	${CC} ${CFLAGS} $(LDFLAGS) -c ${simulate} ${LDFLAGS} -o simulate.o

estim : estim.o optim.o simulate.o model.o $(intermediateC)
	${CC} ${CFLAGS} ${LDFLAGS} $(intermediateC) estim.o optim.o model.o simulate.o ${LDFLAGS} -o estim

monte : optim.o monte.o simulate.o model.o $(intermediateC)
	${MP} ${CFLAGS} ${LDFLAGS} $(intermediateC) optim.o monte.o model.o simulate.o ${LDFLAGS} -o monte

mmic :
	${CC} ${CFLAGS} -c -fopenmp solve_crs_omp.cpp ${LDFLAGS} -o solve_sw_mp_mic.o
	${CC} ${CFLAGS} -c epf.cpp ${LDFLAGS} -o epf_mic.o
	${CC} ${CFLAGS} -c matrixu.cpp -o matrixu_mic.o
	${CC} ${CFLAGS} -c tauchen.cpp -o tachen_mic.o
	${CC} ${CFLAGS} -c utilities.cpp -o utilities_mic.o
	${CC} ${CFLAGS} -c neldmead.cpp -o neldmead_mic.o
	${CC} ${CFLAGS} -c regress.cpp ${LDFLAGS} -o regress_mic.o
	${CC} ${CFLAGS} -c objec.cpp -o objec_mic.o
	${MP} -mmic ${CFLAGS} -c estim_epfq.cpp -o estim_epfq_mic.o

wipe :
	-rm *.o
	-rm *~
	-rm

packup :
	mv -f *.csv Tests
	mv -f *.png Tests
	
