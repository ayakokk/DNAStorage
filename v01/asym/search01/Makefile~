all: bmatrix.o func.o CodeParam.o GenBcode PrintBinCB ConvCheck

bmatrix.o: bmatrix.cpp bmatrix.hpp Makefile
	g++ -Wall -O3 -c -o bmatrix.o bmatrix.cpp

func.o: func.cpp func.hpp Makefile
	g++ -Wall -O3 -c -o func.o func.cpp

CodeParam.o: CodeParam.cpp CodeParam.hpp func.hpp Makefile
	g++ -Wall -O3 -c -o CodeParam.o CodeParam.cpp

GenBcode: GenBcode.cpp bmatrix.o func.o Makefile
	g++ -Wall -O3 -o GenBcode bmatrix.o func.o GenBcode.cpp

PrintBinCB: PrintBinCB.cpp bmatrix.o func.o Makefile
	g++ -Wall -O3 -o PrintBinCB bmatrix.o func.o PrintBinCB.cpp

ConvCheck: ConvCheck.cpp func.o Makefile
	g++ -Wall -O3 -o ConvCheck func.o ConvCheck.cpp

clean:
	rm -f *.o GenBcode PrintBinCB ConvCheck
