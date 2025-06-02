all: bmatrix.o func.o misc.o NBIDSchannel.o InnerCodebook.o Gtable.o FBalg4.o GenCB_ex1 GenCB_ex2a GenCB_ex2b GenCB GenGtable PrintICB PrintGtable SimInner BoundConstSeg

bmatrix.o: bmatrix.cpp bmatrix.hpp Makefile
	g++ -Wall -O3 -c -o bmatrix.o bmatrix.cpp

func.o: func.cpp func.hpp bmatrix.hpp Makefile
	g++ -Wall -O3 -c -o func.o func.cpp

misc.o: misc.cpp misc.hpp Makefile
	g++ -Wall -O3 -c -o misc.o misc.cpp

NBIDSchannel.o: NBIDSchannel.cpp NBIDSchannel.hpp Makefile
	g++ -Wall -O3 -c -o NBIDSchannel.o NBIDSchannel.cpp

InnerCodebook.o: InnerCodebook.cpp InnerCodebook.hpp bmatrix.hpp func.hpp Makefile
	g++ -Wall -O3 -c -o InnerCodebook.o InnerCodebook.cpp

Gtable.o: Gtable.cpp Gtable.hpp bmatrix.hpp func.hpp Makefile
	g++ -Wall -O3 -std=c++11 -c -o Gtable.o Gtable.cpp

FBalg4.o: FBalg4.cpp FBalg4.hpp bmatrix.hpp func.hpp Gtable.hpp Makefile
	g++ -Wall -O3 -std=c++11 -c -o FBalg4.o FBalg4.cpp

GenCB_ex1: GenCB_ex1.cpp bmatrix.o func.o InnerCodebook.o Makefile
	g++ -Wall -O3 -o GenCB_ex1 bmatrix.o func.o InnerCodebook.o GenCB_ex1.cpp

GenCB_ex2a: GenCB_ex2a.cpp bmatrix.o func.o InnerCodebook.o Makefile
	g++ -Wall -O3 -o GenCB_ex2a bmatrix.o func.o InnerCodebook.o GenCB_ex2a.cpp

GenCB_ex2b: GenCB_ex2b.cpp bmatrix.o func.o InnerCodebook.o Makefile
	g++ -Wall -O3 -o GenCB_ex2b bmatrix.o func.o InnerCodebook.o GenCB_ex2b.cpp

GenCB: GenCB.cpp bmatrix.o func.o InnerCodebook.o Makefile
	g++ -Wall -O3 -o GenCB bmatrix.o func.o InnerCodebook.o GenCB.cpp

GenGtable: GenGtable.cpp bmatrix.o func.o Gtable.o Makefile
	g++ -Wall -O3 -o GenGtable bmatrix.o func.o Gtable.o GenGtable.cpp

PrintICB: PrintICB.cpp bmatrix.o func.o InnerCodebook.o Makefile
	g++ -Wall -O3 -o PrintICB bmatrix.o func.o InnerCodebook.o PrintICB.cpp

PrintGtable: PrintGtable.cpp bmatrix.o func.o Gtable.o Makefile
	g++ -Wall -O3 -o PrintGtable bmatrix.o func.o Gtable.o PrintGtable.cpp

SimInner: SimInner.cpp bmatrix.o func.o misc.o NBIDSchannel.o InnerCodebook.o Gtable.o FBalg4.o Makefile
	g++ -Wall -O3 -o SimInner bmatrix.o func.o misc.o NBIDSchannel.o InnerCodebook.o Gtable.o FBalg4.o SimInner.cpp

BoundConstSeg: BoundConstSeg.cpp bmatrix.o func.o misc.o Makefile
	g++ -Wall -O3 -o BoundConstSeg bmatrix.o func.o misc.o BoundConstSeg.cpp

clean:
	rm -f *.o GenCB_ex1 GenCB_ex2a GenCB_ex2b GenCB GenGtable PrintICB PrintGtable SimInner BoundConstSeg
