all: func.o MutualInfo.o InnerCodebook.o ReadICB simEnc

func.o: func.cpp func.hpp Makefile
	g++ -Wall -O3 -c -o func.o func.cpp

MutualInfo.o: MutualInfo.cpp MutualInfo.hpp Makefile
	g++ -Wall -O3 -c -o MutualInfo.o MutualInfo.cpp

InnerCodebook.o: InnerCodebook.cpp InnerCodebook.hpp Makefile
	g++ -Wall -O3 -c -o InnerCodebook.o InnerCodebook.cpp

ReadICB: ReadICB.cpp InnerCodebook.o Makefile
	g++ -Wall -O3 -o ReadICB InnerCodebook.o ReadICB.cpp

simEnc: simEnc.cpp func.o MutualInfo.o InnerCodebook.o Makefile
	g++ -Wall -O3 -o simEnc func.o MutualInfo.o InnerCodebook.o simEnc.cpp

clean:
	rm -f ReadICB simEnc *.o

