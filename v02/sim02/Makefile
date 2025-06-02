CXX = g++
CXXFLAGS = -Wall -O3 -I/Users/ayako/.julia/juliaup/julia-1.11.5+0.aarch64.apple.darwin14/include/julia
LDFLAGS = -L/Users/ayako/.julia/juliaup/julia-1.11.5+0.aarch64.apple.darwin14/lib -ljulia

all: func.o ChannelMatrix.o InnerCodebook.o IDSchannel.o SLFBAdec.o ATGCConverter.o JuliaChannel.o CBsearch GenEncCM sim sim-new-channel

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

# 新規追加
ATGCConverter.o: ATGCConverter.cpp ATGCConverter.hpp Makefile
	$(CXX) $(CXXFLAGS) -c -o ATGCConverter.o ATGCConverter.cpp

# simターゲットを修正
sim: sim.cpp func.o ChannelMatrix.o InnerCodebook.o IDSchannel.o SLFBAdec.o ATGCConverter.o JuliaChannel.o Makefile
	$(CXX) $(CXXFLAGS) -o sim func.o ChannelMatrix.o InnerCodebook.o IDSchannel.o SLFBAdec.o sim.cpp $(LDFLAGS)

GenEncCM: GenEncCM.cpp func.o ChannelMatrix.o InnerCodebook.o Makefile
	$(CXX) $(CXXFLAGS) -o GenEncCM func.o ChannelMatrix.o InnerCodebook.o GenEncCM.cpp $(LDFLAGS)

CBsearch: CBsearch.cpp func.o Makefile
	$(CXX) $(CXXFLAGS) -o CBsearch func.o CBsearch.cpp $(LDFLAGS)

# 新規追加  
JuliaChannel.o: JuliaChannel.cpp JuliaChannel.hpp ATGCConverter.hpp Makefile
	$(CXX) $(CXXFLAGS) -c -o JuliaChannel.o JuliaChannel.cpp

sim-new-channel: sim-new-channel.cpp func.o ChannelMatrix.o InnerCodebook.o IDSchannel.o SLFBAdec.o ATGCConverter.o JuliaChannel.o Makefile
	$(CXX) $(CXXFLAGS) -o sim-new-channel func.o ChannelMatrix.o InnerCodebook.o IDSchannel.o SLFBAdec.o ATGCConverter.o JuliaChannel.o sim-new-channel.cpp $(LDFLAGS)




normal: sim
	@echo "Starting normal experiments..."
	@mkdir -p log/Normal
	sh/Normal/ex01a.sh
	sh/Normal/ex01b.sh
	sh/Normal/ex01c.sh
	@echo "All normal experiments completed!"

# 並列実行版（バックグラウンド実行）
extreme-parallel: sim
	@echo "Starting parallel Extreme experiments..."
	@mkdir -p log/Extreme
	sh/Extreme/ex01a.sh & \
	sh/Extreme/ex01b.sh & \
	sh/Extreme/ex01c.sh & \
	wait
	@echo "All parallel Extreme experiments completed!"

# Extreme_reverse実験を並列実行
extreme-reverse-parallel: sim
	@echo "Starting parallel Extreme_reverse experiments..."
	@mkdir -p log
	sh/Extreme_reverse/ex01a.sh & \
	sh/Extreme_reverse/ex01b.sh & \
	sh/Extreme_reverse/ex01c.sh & \
	wait
	@echo "All parallel Extreme_reverse experiments completed!"

light-parallel : sim
	@echo "Starting light experiments..."
	@mkdir -p log/Light
	sh/Light/ex01a.sh & \
	sh/Light/ex01b.sh & \
	sh/Light/ex01c.sh & \
	wait
	@echo "All light experiments completed!"

light-reverse-parallel : sim
	@echo "Starting light_reverse experiments..."
	@mkdir -p log/Light_reverse
	sh/Light_reverse/ex01a.sh & \
	sh/Light_reverse/ex01b.sh & \
	sh/Light_reverse/ex01c.sh & \
	wait
	@echo "All light_reverse experiments completed!"

medium-parallel : sim
	@echo "Starting medium experiments..."
	@mkdir -p log/Medium
	sh/Medium/ex01a.sh & \
	sh/Medium/ex01b.sh & \
	sh/Medium/ex01c.sh & \
	wait
	@echo "All medium experiments completed!"

medium-reverse-parallel : sim
	@echo "Starting medium_reverse experiments..."
	@mkdir -p log/Medium_reverse
	sh/Medium_reverse/ex01a.sh & \
	sh/Medium_reverse/ex01b.sh & \
	sh/Medium_reverse/ex01c.sh & \
	wait
	@echo "All medium_reverse experiments completed!"

strong-parallel : sim
	@echo "Starting strong experiments..."
	@mkdir -p log/Strong
	sh/Strong/ex01a.sh & \
	sh/Strong/ex01b.sh & \
	sh/Strong/ex01c.sh & \
	wait
	@echo "All strong experiments completed!"

strong-reverse-parallel : sim
	@echo "Starting strong_reverse experiments..."
	@mkdir -p log/Strong_reverse
	sh/Strong_reverse/ex01a.sh & \
	sh/Strong_reverse/ex01b.sh & \
	sh/Strong_reverse/ex01c.sh & \
	wait
	@echo "All strong_reverse experiments completed!"


clean:
	rm -f CBsearch GenEncCM sim *.o

