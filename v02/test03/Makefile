all: InnerCodebook.o ReadICB

InnerCodebook.o: InnerCodebook.cpp InnerCodebook.hpp Makefile
	g++ -Wall -O3 -c -o InnerCodebook.o InnerCodebook.cpp

ReadICB: ReadICB.cpp InnerCodebook.o Makefile
	g++ -Wall -O3 -o ReadICB InnerCodebook.o ReadICB.cpp

clean:
	rm -f ReadICB *.o

