# C++コンパイラを指定
CXX = g++

# デバッグ版（開発用・OpenMPなし）
CXXFLAGS_DEBUG = -Wall -g -fsanitize=address -std=c++11
LDFLAGS_DEBUG =

# OpenMP版（デバッグ）
CXXFLAGS_DEBUG_OMP = -Wall -g -fsanitize=address -std=c++11 -fopenmp
LDFLAGS_DEBUG_OMP = -fopenmp

# 高性能版（TSUBAME用）
CXXFLAGS_RELEASE = -Wall -O3 -std=c++11 -fopenmp -DNDEBUG -march=native
LDFLAGS_RELEASE = -fopenmp

# デフォルトはデバッグ版（OpenMPなし）
CXXFLAGS = $(CXXFLAGS_DEBUG)
LDFLAGS = $(LDFLAGS_DEBUG)

# 最終的に生成する実行ファイルを指定
TARGETS = CBsearch GenEncCM sim test_tables

# 'make' とだけ打った時に 'all' が実行される
all: $(TARGETS)

# OpenMP版ビルド（ローカルテスト用）
openmp: clean
	$(MAKE) CXXFLAGS="$(CXXFLAGS_DEBUG_OMP)" LDFLAGS="$(LDFLAGS_DEBUG_OMP)" $(TARGETS)
	@echo "# OpenMP build completed for local testing"

# TSUBAME用高性能版ビルド
release: clean
	$(MAKE) CXXFLAGS="$(CXXFLAGS_RELEASE)" LDFLAGS="$(LDFLAGS_RELEASE)" $(TARGETS)
	@echo "# High-performance build completed for TSUBAME"

# OpenMPテスト用（小規模）
test-openmp: openmp
	@echo "# Testing OpenMP implementation..."
	OMP_NUM_THREADS=4 ./sim ICB/ex01-new-4gen/ 3 0

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
	$(CXX) $(CXXFLAGS) -o CBsearch func.o CBsearch.cpp $(LDFLAGS)

GenEncCM: GenEncCM.cpp func.o ChannelMatrix.o InnerCodebook.o Makefile
	$(CXX) $(CXXFLAGS) -o GenEncCM func.o ChannelMatrix.o InnerCodebook.o GenEncCM.cpp $(LDFLAGS)

sim: sim.cpp $(OBJS) Makefile
	$(CXX) $(CXXFLAGS) -o sim $(OBJS) sim.cpp $(LDFLAGS)

test_tables: test_tables.cpp $(OBJS) Makefile
	$(CXX) $(CXXFLAGS) -o test_tables $(OBJS) test_tables.cpp $(LDFLAGS)

# クリーンアップルール
clean:
	rm -f $(TARGETS) *.o