all: func.o ChannelMatrix.o InnerCodebook.o IDSchannel.o SLFBAdec.o CBsearch GenEncCM sim

func.o: func.cpp func.hpp Makefile
	g++ -Wall -O3 -c -o func.o func.cpp

ChannelMatrix.o: ChannelMatrix.cpp ChannelMatrix.hpp Makefile
	g++ -Wall -O3 -c -o ChannelMatrix.o ChannelMatrix.cpp

InnerCodebook.o: InnerCodebook.cpp InnerCodebook.hpp Makefile
	g++ -Wall -O3 -c -o InnerCodebook.o InnerCodebook.cpp

IDSchannel.o: IDSchannel.cpp IDSchannel.hpp Makefile
	g++ -Wall -O3 -c -o IDSchannel.o IDSchannel.cpp

SLFBAdec.o: SLFBAdec.cpp SLFBAdec.hpp InnerCodebook.hpp ChannelMatrix.hpp IDSchannel.hpp Makefile
	g++ -Wall -O3 -c -o SLFBAdec.o SLFBAdec.cpp

CBsearch: CBsearch.cpp func.o Makefile
	g++ -Wall -O3 -o CBsearch func.o CBsearch.cpp

GenEncCM: GenEncCM.cpp func.o ChannelMatrix.o InnerCodebook.o Makefile
	g++ -Wall -O3 -o GenEncCM func.o ChannelMatrix.o InnerCodebook.o GenEncCM.cpp

sim: sim.cpp func.o ChannelMatrix.o InnerCodebook.o IDSchannel.o SLFBAdec.o Makefile
	g++ -Wall -O3 -o sim func.o ChannelMatrix.o InnerCodebook.o IDSchannel.o SLFBAdec.o sim.cpp

clean:
	rm -f CBsearch GenEncCM sim *.o

