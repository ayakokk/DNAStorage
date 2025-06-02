all: bmatrix.o func.o CodeParam.o GenBCB RunLenBCB ConvCheck

bmatrix.o: bmatrix.cpp bmatrix.hpp Makefile
	g++ -Wall -O3 -c -o bmatrix.o bmatrix.cpp

func.o: func.cpp func.hpp Makefile
	g++ -Wall -O3 -c -o func.o func.cpp

CodeParam.o: CodeParam.cpp CodeParam.hpp func.hpp Makefile
	g++ -Wall -O3 -c -o CodeParam.o CodeParam.cpp

GenBCB: GenBCB.cpp bmatrix.o func.o CodeParam.o Makefile
	g++ -Wall -O3 -o GenBCB bmatrix.o func.o CodeParam.o GenBCB.cpp

RunLenBCB: RunLenBCB.cpp bmatrix.o func.o Makefile
	g++ -Wall -O3 -o RunLenBCB bmatrix.o func.o RunLenBCB.cpp

ConvCheck: ConvCheck.cpp func.o Makefile
	g++ -Wall -O3 -o ConvCheck func.o ConvCheck.cpp

clean:
	rm -f *.o GenBCB PrintBinCB RunLenBCB ConvCheck
