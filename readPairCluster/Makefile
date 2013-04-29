CFLAGS=-I /scr/talkowski/tools/bin/samtools-current -L /scr/talkowski/tools/bin/samtools-current -lbam -L ~/local/lib -lz

readpairCluster: cluster.o readPair.o
	g++ -g -o readPairCluster cluster.o readPair.o $(CFLAGS)

cluster.o: cluster.cpp simpleGraph.hpp
	g++ -g -c cluster.cpp $(CFLAGS)

readPair.o: readPair.cpp readPair.hpp 
	g++ -g -c readPair.cpp $(CFLAGS)
