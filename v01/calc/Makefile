all: RunLenRnd GCB_Rnd

RunLenRnd: RunLenRnd.cpp Makefile
	g++ -Wall -O3 -o RunLenRnd RunLenRnd.cpp

GCB_Rnd: GCB_Rnd.cpp Makefile
	g++ -Wall -O3 -o GCB_Rnd GCB_Rnd.cpp

clean:
	rm -f *.o RunLenRnd GCB_Rnd
