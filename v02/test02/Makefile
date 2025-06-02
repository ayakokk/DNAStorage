all: func.o InnerCodec.o sim

func.o: func.cpp func.hpp Makefile
	g++ -Wall -O3 -c -o func.o func.cpp

InnerCodec.o: InnerCodec.cpp InnerCodecFenc.cpp InnerCodec.hpp Makefile
	g++ -Wall -O3 -c -o InnerCodec.o InnerCodec.cpp

sim: sim.cpp func.o InnerCodec.o Makefile
	g++ -Wall -O3 -o sim func.o InnerCodec.o sim.cpp

clean:
	rm -f *.o sim
