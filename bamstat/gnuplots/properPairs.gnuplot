set style histogram rowstacked
set boxwidth 0.9 relative
set style fill solid 0.75
set title "Number of Mapped Read-Pairs by Orientation"
set xlabel "Read Pair Orientation"
set ylabel "Number of Read Pairs"
set key left
plot "test.dat" using 2:xticlabels(1) title "Proper Pairs" lt rgb "sea-green", \
 "test.dat" using 3:xticlabels(1) lt rgb "gray50" title "Not Properly Paired"
