CC = g++-mp-4.5
EXE = myFEM.exe
SRC = myFEM.cpp
OBJ = myFEM.o

LDLIBS = -I/Users/dalg24/Packages/SuiteSparse/UFconfig \
	       -I/Users/dalg24/Packages/SuiteSparse/UMFPACK/Include \
		     -I/Users/dalg24/Packages/SuiteSparse/AMD/Include \
         -I/Users/dalg24/myFEM/include
LDFLAGS = -L/opt/local/lib \
	        -L/Users/dalg24/Packages/lapack/lapack-3.3.1 \
	        -L/Users/dalg24/Packages/SuiteSparse/UMFPACK/Lib \
	        -L/Users/dalg24/Packages/SuiteSparse/AMD/Lib \
          -lm -lamd -lumfpack -lmetis -lblas
CFLAGS =  -g -Wall
CDEFS = -DEBILE

$(EXE) : $(OBJ)
	$(CC) $(OBJ) $(LDFLAGS) -o $(EXE)

$(OBJ) : $(SRC) 
	$(CC) $(CDEFS) ${CFLAGS} $(LDLIBS) -c $(SRC)

build-tests:
	g++ testQuadratureRule.cpp -o testQuadratureRule.exe -Wall -I/Users/dalg24/myFEM/include
	g++ testBasisFunctions.cpp -o testBasisFunctions.exe -Wall -I/Users/dalg24/myFEM/include

clean:
	rm $(OBJ)

cleanall:
	rm $(OBJ)
	rm -f $(EXE)

