# C++コンパイラを指定
CXX = g++
# コンパイルフラグを指定
# メモ化アプローチにより大幅なメモリ効率化を実現（40.5GB → 数KB）
CXXFLAGS = -Wall -g -fsanitize=address -std=c++11

# 最終的に生成する実行ファイルを指定
TARGETS = CBsearch GenEncCM sim test_tables

# 'make' とだけ打った時に 'all' が実行される
all: $(TARGETS)

# オブジェクトファイルの依存関係
OBJS = func.o ChannelMatrix.o InnerCodebook.o IDSchannel.o SLFBAdec.o

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

# 実行ファイルの生成ルール
CBsearch: CBsearch.cpp func.o Makefile
	$(CXX) $(CXXFLAGS) -o CBsearch func.o CBsearch.cpp

GenEncCM: GenEncCM.cpp func.o ChannelMatrix.o InnerCodebook.o Makefile
	$(CXX) $(CXXFLAGS) -o GenEncCM func.o ChannelMatrix.o InnerCodebook.o GenEncCM.cpp

sim: sim.cpp $(OBJS) Makefile
	$(CXX) $(CXXFLAGS) -o sim $(OBJS) sim.cpp

test_tables: test_tables.cpp $(OBJS) Makefile
	$(CXX) $(CXXFLAGS) -o test_tables $(OBJS) test_tables.cpp

# クリーンアップルール
clean:
	rm -f $(TARGETS) *.o