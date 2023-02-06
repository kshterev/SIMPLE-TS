## Compilers
#CC = icpc
#CC = g++
CC = mpic++

SIMPLE_TS_MPI: SIMPLE_TS_MPI.cpp Circles_Objects.cpp Polyhedrons_Objects.cpp Objects_for_Solving.cpp DomainDecomposition2D.cpp
## Options for intel compiler (icpc):
#	$(CC) -xSSE2 -O3 -static -ipo SIMPLE_TS_MPI.cpp Circles_Objects.cpp Polyhedrons_Objects.cpp Objects_for_Solving.cpp DomainDecomposition2D.cpp -o SIMPLE_TS_without_MPI_static
## Options for gcc compiler (g++):
#	$(CC) -msse2 -O3 -static SIMPLE_TS_MPI.cpp Circles_Objects.cpp Polyhedrons_Objects.cpp Objects_for_Solving.cpp DomainDecomposition2D.cpp -o SIMPLE_TS_without_MPI_static
## Options for openmpi compiler (mpic++):
	$(CC) -msse2 -O3 SIMPLE_TS_MPI.cpp Circles_Objects.cpp Polyhedrons_Objects.cpp Objects_for_Solving.cpp DomainDecomposition2D.cpp -o SIMPLE_TS_MPI



## Delete *.o files:
	rm -f Circles_Objects.o
	rm -f DomainDecomposition2D.o
	rm -f Objects_for_Solving.o
	rm -f Polyhedrons_Objects.o
	rm -f SIMPLE_TS_MPI.o
