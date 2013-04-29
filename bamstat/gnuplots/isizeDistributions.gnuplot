
set boxwidth 0.9 relative
set style fill solid 0.75
set title 
plot "longReads.histdata" using 2:xticlabels(3) title "Counts (Proper R-F pairs)" with boxes lt rgb "sea-green"
