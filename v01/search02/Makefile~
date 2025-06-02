all: func.o func2.o BaseFile.o combin.o EnumVect.o bmatrix.o NBIDSchannel.o InnerCodebook.o Gtable.o FBalg4.o StatICB.o GenBinCode ComBinCode GenICB2 SimInner 

func.o: func.cpp func.hpp Makefile
	g++ -Wall -O3 -c -o func.o func.cpp

func2.o: func2.cpp func2.hpp Makefile
	g++ -Wall -O3 -c -o func2.o func2.cpp

BaseFile.o: BaseFile.cpp BaseFile.hpp Makefile
	g++ -Wall -O3 -std=c++14 -c -o BaseFile.o BaseFile.cpp

combin.o: combin.cpp combin.hpp Makefile
	g++ -Wall -O3 -std=c++14 -c -o combin.o combin.cpp

EnumVect.o: EnumVect.cpp EnumVect.hpp Makefile
	g++ -Wall -O3 -std=c++14 -c -o EnumVect.o EnumVect.cpp

bmatrix.o: bmatrix.cpp bmatrix.hpp Makefile
	g++ -Wall -O3 -c -o bmatrix.o bmatrix.cpp

NBIDSchannel.o: NBIDSchannel.cpp NBIDSchannel.hpp Makefile
	g++ -Wall -O3 -c -o NBIDSchannel.o NBIDSchannel.cpp

InnerCodebook.o: InnerCodebook.cpp InnerCodebook.hpp bmatrix.hpp func.hpp Makefile
	g++ -Wall -O3 -c -o InnerCodebook.o InnerCodebook.cpp

Gtable.o: Gtable.cpp Gtable.hpp bmatrix.hpp func.hpp Makefile
	g++ -Wall -O3 -std=c++14 -c -o Gtable.o Gtable.cpp

FBalg4.o: FBalg4.cpp FBalg4.hpp bmatrix.hpp func.hpp Gtable.hpp Makefile
	g++ -Wall -O3 -std=c++14 -c -o FBalg4.o FBalg4.cpp

StatICB.o: StatICB.cpp StatICB.hpp func.hpp Makefile
	g++ -Wall -O3 -std=c++14 -c -o StatICB.o StatICB.cpp

GenBinCode: GenBinCode.cpp func2.o combin.o bmatrix.o Makefile
	g++ -Wall -O3 -std=c++14 -o GenBinCode func2.o combin.o bmatrix.o GenBinCode.cpp

ComBinCode: ComBinCode.cpp func2.o BaseFile.o EnumVect.o Makefile
	g++ -Wall -O3 -std=c++14 -o ComBinCode func2.o BaseFile.o EnumVect.o ComBinCode.cpp

GenICB2: GenICB2.cpp bmatrix.o func.o InnerCodebook.o Makefile
	g++ -Wall -O3 -std=c++14 -o GenICB2 bmatrix.o func.o InnerCodebook.o GenICB2.cpp

SimInner: SimInner.cpp bmatrix.o func.o NBIDSchannel.o InnerCodebook.o Gtable.o FBalg4.o StatICB.o Makefile
	g++ -Wall -O3 -o SimInner bmatrix.o func.o NBIDSchannel.o InnerCodebook.o Gtable.o FBalg4.o StatICB.o SimInner.cpp

clean:
	rm -f *.o GenBinCode ComBinCode GenICB2 SimInner

