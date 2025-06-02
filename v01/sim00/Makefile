all: bmatrix.o func.o NBIDSchannel.o InnerCodebook.o Gtable.o FBalg4.o GenICB GenICBmk GenGtable PrintICB PrintGtable SimInner CalcIxy BoundConstSeg RL_LGCB

bmatrix.o: bmatrix.cpp bmatrix.hpp Makefile
	g++ -Wall -O3 -c -o bmatrix.o bmatrix.cpp

func.o: func.cpp func.hpp bmatrix.hpp Makefile
	g++ -Wall -O3 -c -o func.o func.cpp

NBIDSchannel.o: NBIDSchannel.cpp NBIDSchannel.hpp Makefile
	g++ -Wall -O3 -c -o NBIDSchannel.o NBIDSchannel.cpp

InnerCodebook.o: InnerCodebook.cpp InnerCodebook.hpp bmatrix.hpp func.hpp Makefile
	g++ -Wall -O3 -c -o InnerCodebook.o InnerCodebook.cpp

Gtable.o: Gtable.cpp Gtable.hpp bmatrix.hpp func.hpp Makefile
	g++ -Wall -O3 -std=c++11 -c -o Gtable.o Gtable.cpp

FBalg4.o: FBalg4.cpp FBalg4.hpp bmatrix.hpp func.hpp Gtable.hpp Makefile
	g++ -Wall -O3 -std=c++11 -c -o FBalg4.o FBalg4.cpp

GenICB: GenICB.cpp bmatrix.o func.o InnerCodebook.o Makefile
	g++ -Wall -O3 -o GenICB bmatrix.o func.o InnerCodebook.o GenICB.cpp

GenICBmk: GenICBmk.cpp bmatrix.o func.o InnerCodebook.o Makefile
	g++ -Wall -O3 -o GenICBmk bmatrix.o func.o InnerCodebook.o GenICBmk.cpp

GenGtable: GenGtable.cpp bmatrix.o func.o Gtable.o Makefile
	g++ -Wall -O3 -o GenGtable bmatrix.o func.o Gtable.o GenGtable.cpp

PrintICB: PrintICB.cpp bmatrix.o func.o InnerCodebook.o Makefile
	g++ -Wall -O3 -o PrintICB bmatrix.o func.o InnerCodebook.o PrintICB.cpp

PrintGtable: PrintGtable.cpp bmatrix.o func.o Gtable.o Makefile
	g++ -Wall -O3 -o PrintGtable bmatrix.o func.o Gtable.o PrintGtable.cpp

SimInner: SimInner.cpp bmatrix.o func.o NBIDSchannel.o InnerCodebook.o Gtable.o FBalg4.o Makefile
	g++ -Wall -O3 -o SimInner bmatrix.o func.o NBIDSchannel.o InnerCodebook.o Gtable.o FBalg4.o SimInner.cpp

CalcIxy: CalcIxy.cpp bmatrix.o func.o NBIDSchannel.o InnerCodebook.o Gtable.o FBalg4.o Makefile
	g++ -Wall -O3 -o CalcIxy bmatrix.o func.o NBIDSchannel.o InnerCodebook.o Gtable.o FBalg4.o CalcIxy.cpp

BoundConstSeg: BoundConstSeg.cpp bmatrix.o func.o Makefile
	g++ -Wall -O3 -o BoundConstSeg bmatrix.o func.o BoundConstSeg.cpp

RL_LGCB: RL_LGCB.cpp bmatrix.o func.o InnerCodebook.o Makefile
	g++ -Wall -O3 -o RL_LGCB bmatrix.o func.o InnerCodebook.o RL_LGCB.cpp

clean:
	rm -f *.o GenICB GenICBmk GenGtable PrintICB PrintGtable SimInner CalcIxy BoundConstSeg RL_LGCB
