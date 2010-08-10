#

L = -L/usr/lib
I = -I/usr/include/suitesparse
CFLAGS =  -lm -lamd -lumfpack

#CDEFS = -DONED
CDEFS = -DTWOD

all:
	g++ -c -O3 $(CDEFS) $(CFLAGS) $(I) coomatrix.cpp myFEsolver.cpp
	g++ -o avanti $(I) $(L) $(CFLAGS) coomatrix.o myFEsolver.o

clean:
	rm *.o
	rm *~

cleanall:
	rm *.o
	rm *~
	rm -f avanti


