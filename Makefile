CC = g++-mp-4.5
SRC = myFEM.cpp
EXE = myFEM.exe

I = -I/Users/dalg24/Packages/SuiteSparse/UFconfig \
	  -I/Users/dalg24/Packages/SuiteSparse/UMFPACK/Include \
		-I/Users/dalg24/Packages/SuiteSparse/AMD/Include

L = -L/opt/local/lib \
		-L/Users/dalg24/Packages/lapack/lapack-3.3.1 \
		-L/Users/dalg24/Packages/SuiteSparse/UMFPACK/Lib \
		-L/Users/dalg24/Packages/SuiteSparse/AMD/Lib
CFLAGS = -lm -lamd -lumfpack -lmetis -lblas

CDEFS = -DEBILE
LDFLAGS = -arch x86_64

all:
	${CC} $(CDEFS) $(I) -c ${SRC} 
	${CC} *.o $(L) $(CFLAGS) $(LDFLAGS) -o ${EXE}

clean:
	rm *.o

cleanall:
	rm *.o
	rm -f ${EXE}

