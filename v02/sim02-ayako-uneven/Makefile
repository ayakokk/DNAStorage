CXX = g++
CXXFLAGS = -Wall -O3

all: func.o ChannelMatrix.o InnerCodebook.o IDSchannel.o SLFBAdec.o ATGCConverter.o CBsearch GenEncCM sim sim-new-channel

func.o: func.cpp func.hpp Makefile
	$(CXX) $(CXXFLAGS) -c -o func.o func.cpp

ChannelMatrix.o: ChannelMatrix.cpp ChannelMatrix.hpp Makefile
	$(CXX) $(CXXFLAGS) -c -o ChannelMatrix.o ChannelMatrix.cpp

InnerCodebook.o: InnerCodebook.cpp InnerCodebook.hpp Makefile
	$(CXX) $(CXXFLAGS) -c -o InnerCodebook.o InnerCodebook.cpp

IDSchannel.o: IDSchannel.cpp IDSchannel.hpp Makefile
	$(CXX) $(CXXFLAGS) -c -o IDSchannel.o IDSchannel.cpp

SLFBAdec.o: SLFBAdec.cpp SLFBAdec.hpp InnerCodebook.hpp ChannelMatrix.hpp IDSchannel.hpp Makefile
	$(CXX) $(CXXFLAGS) -c -o SLFBAdec.o SLFBAdec.cpp

ATGCConverter.o: ATGCConverter.cpp ATGCConverter.hpp Makefile
	$(CXX) $(CXXFLAGS) -c -o ATGCConverter.o ATGCConverter.cpp

# 元のsim
sim: sim.cpp func.o ChannelMatrix.o InnerCodebook.o IDSchannel.o SLFBAdec.o ATGCConverter.o Makefile
	$(CXX) $(CXXFLAGS) -o sim func.o ChannelMatrix.o InnerCodebook.o IDSchannel.o SLFBAdec.o ATGCConverter.o sim.cpp

# 新しいsim-new-channel
sim-new-channel: sim-new-channel.cpp func.o ChannelMatrix.o InnerCodebook.o IDSchannel.o SLFBAdec.o ATGCConverter.o Makefile
	$(CXX) $(CXXFLAGS) -o sim-new-channel func.o ChannelMatrix.o InnerCodebook.o IDSchannel.o SLFBAdec.o ATGCConverter.o sim-new-channel.cpp

GenEncCM: GenEncCM.cpp func.o ChannelMatrix.o InnerCodebook.o Makefile
	$(CXX) $(CXXFLAGS) -o GenEncCM func.o ChannelMatrix.o InnerCodebook.o GenEncCM.cpp

CBsearch: CBsearch.cpp func.o Makefile
	$(CXX) $(CXXFLAGS) -o CBsearch func.o CBsearch.cpp

clean:
	rm -f CBsearch GenEncCM sim sim-new-channel *.o