all: func.o combin.o MetricCalc.o bmatrix.o GenBinCode RndBinCode

func.o: func.cpp func.hpp Makefile
	g++ -Wall -O3 -c -o func.o func.cpp

combin.o: combin.cpp combin.hpp Makefile
	g++ -Wall -O3 -std=c++14 -c -o combin.o combin.cpp

MetricCalc.o: MetricCalc.cpp MetricCalc.hpp Makefile
	g++ -Wall -O3 -c -o MetricCalc.o MetricCalc.cpp

bmatrix.o: bmatrix.cpp bmatrix.hpp Makefile
	g++ -Wall -O3 -c -o bmatrix.o bmatrix.cpp

GenBinCode: GenBinCode.cpp func.o combin.o MetricCalc.o bmatrix.o Makefile
	g++ -Wall -O3 -std=c++14 -o GenBinCode func.o combin.o MetricCalc.o bmatrix.o GenBinCode.cpp

RndBinCode: RndBinCode.cpp func.o combin.o MetricCalc.o bmatrix.o Makefile
	g++ -Wall -O3 -std=c++14 -o RndBinCode func.o combin.o MetricCalc.o bmatrix.o RndBinCode.cpp

clean:
	rm -f *.o GenBinCode RndBinCode

