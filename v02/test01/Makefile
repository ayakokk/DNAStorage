all: InnerCodec.o genIcodec

InnerCodec.o: InnerCodec.cpp InnerCodec.hpp Makefile
	g++ -Wall -O3 -c -o InnerCodec.o InnerCodec.cpp

genIcodec: genIcodec.cpp InnerCodec.o Makefile
	g++ -Wall -O3 -o genIcodec InnerCodec.o genIcodec.cpp

clean:
	rm -f *.o genIcodec
