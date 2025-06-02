all: func.o bmatrix.o InnerCodebook.o SearchICB

func.o: func.cpp func.hpp Makefile
	g++ -Wall -O3 -std=c++14 -c -o func.o func.cpp

bmatrix.o: bmatrix.cpp bmatrix.hpp Makefile
	g++ -Wall -O3 -std=c++14 -c -o bmatrix.o bmatrix.cpp

InnerCodebook.o: InnerCodebook.cpp InnerCodebook.hpp bmatrix.hpp func.hpp Makefile
	g++ -Wall -O3 -std=c++14 -c -o InnerCodebook.o InnerCodebook.cpp

SearchICB: SearchICB.cpp bmatrix.o func.o InnerCodebook.o Makefile
	g++ -Wall -O3 -std=c++14 -o SearchICB bmatrix.o func.o InnerCodebook.o SearchICB.cpp

clean:
	rm -f *.o SearchICB

