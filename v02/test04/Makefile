all: func.o InnerCodebook.o ReadICB sim

func.o: func.cpp func.hpp Makefile
	g++ -Wall -O3 -c -o func.o func.cpp

InnerCodebook.o: InnerCodebook.cpp InnerCodebook.hpp Makefile
	g++ -Wall -O3 -c -o InnerCodebook.o InnerCodebook.cpp

ReadICB: ReadICB.cpp InnerCodebook.o Makefile
	g++ -Wall -O3 -o ReadICB InnerCodebook.o ReadICB.cpp

sim: sim.cpp func.o InnerCodebook.o Makefile
	g++ -Wall -O3 -o sim func.o InnerCodebook.o sim.cpp

clean:
	rm -f ReadICB sim *.o

