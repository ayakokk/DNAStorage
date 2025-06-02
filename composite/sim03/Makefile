all: func.o ChannelMatrix.o InnerCodebook.o CIDSchannel.o CSFBAdec.o CPScodec.o ConcatAB.o CBsearch GenEncCM sim

func.o: func.cpp func.hpp Makefile
	g++ -Wall -O3 -c -o func.o func.cpp

ChannelMatrix.o: ChannelMatrix.cpp ChannelMatrix.hpp Makefile
	g++ -Wall -O3 -c -o ChannelMatrix.o ChannelMatrix.cpp

InnerCodebook.o: InnerCodebook.cpp InnerCodebook.hpp Makefile
	g++ -Wall -O3 -c -o InnerCodebook.o InnerCodebook.cpp

CIDSchannel.o: CIDSchannel.cpp CIDSchannel.hpp Makefile
	g++ -Wall -O3 -c -o CIDSchannel.o CIDSchannel.cpp

CSFBAdec.o: CSFBAdec.cpp CSFBAdec.hpp InnerCodebook.hpp ChannelMatrix.hpp Makefile
	g++ -Wall -O3 -c -o CSFBAdec.o CSFBAdec.cpp

CPScodec.o: CPScodec.cpp CPScodec.hpp Makefile
	g++ -Wall -O3 -c -o CPScodec.o CPScodec.cpp

ConcatAB.o: ConcatAB.cpp ConcatAB.hpp InnerCodebook.hpp CSFBAdec.hpp CPScodec.hpp Makefile
	g++ -Wall -O3 -c -o ConcatAB.o ConcatAB.cpp

CBsearch: CBsearch.cpp func.o Makefile
	g++ -Wall -O3 -o CBsearch func.o CBsearch.cpp

GenEncCM: GenEncCM.cpp func.o ChannelMatrix.o InnerCodebook.o Makefile
	g++ -Wall -O3 -o GenEncCM func.o ChannelMatrix.o InnerCodebook.o GenEncCM.cpp

sim: sim.cpp func.o ChannelMatrix.o InnerCodebook.o CIDSchannel.o CSFBAdec.o CPScodec.o ConcatAB.o Makefile
	g++ -Wall -O3 -o sim func.o ChannelMatrix.o InnerCodebook.o CIDSchannel.o CSFBAdec.o CPScodec.o ConcatAB.o sim.cpp

clean:
	rm -f CBsearch GenEncCM sim *.o

