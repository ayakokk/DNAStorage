all: func.o ChannelMatrix.o InnerCodebook.o IDSchannel.o SLFBAdec.o CBsearch GenEncCM sim test test_pyx test_v2

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

test: test.cpp func.o ChannelMatrix.o InnerCodebook.o IDSchannel.o SLFBAdec.o Makefile
	g++ -Wall -O3 -o test func.o ChannelMatrix.o InnerCodebook.o IDSchannel.o SLFBAdec.o test.cpp

test_pyx: test_pyx.cpp func.o ChannelMatrix.o InnerCodebook.o IDSchannel.o SLFBAdec.o Makefile
	g++ -Wall -O3 -o test_pyx func.o ChannelMatrix.o InnerCodebook.o IDSchannel.o SLFBAdec.o test_pyx.cpp

test_v2: test_v2.cpp func.o ChannelMatrix.o InnerCodebook.o IDSchannel.o SLFBAdec.o Makefile
	g++ -Wall -O3 -o test_v2 func.o ChannelMatrix.o InnerCodebook.o IDSchannel.o SLFBAdec.o test_v2.cpp

clean:
	rm -f CBsearch GenEncCM sim test test_pyx test_v2 *.o

