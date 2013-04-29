CFLAGS = -I /scr/talkowski/tools/bin/samtools -L /scr/talkowski/tools/bin/samtools -lbam -L ~/local/lib -lz

bamstat: bamstat.cpp histogram.o 
	g++ -g -o bamstat bamstat.cpp histogram.o discbins.o $(CFLAGS) 

histogram.o: histogram.cpp
	g++ -c histogram.cpp $(CFLAGS)

discbins.o: discbins.cpp
	g++ -c discbins.cpp $(CFLAGS)

bam_select: bam_select.cpp
	g++ -o bam_select bam_select.cpp $(CFLAGS)

