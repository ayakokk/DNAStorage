reset
set grid
set yrange [0.5:3.5]
set key outside right Left
set format y "%.1f"
set xlabel "r"
set ylabel "mutual information"
plot "minfoA.dat" using 1:2  with linespoints lc 1 pt 1 dt 3 title "(6,3) A",\
     "minfoA.dat" using 1:3  with linespoints lc 1 pt 1 dt 2 title "(6,3) B",\
     "minfoA.dat" using 1:4  with linespoints lc 1 pt 1 dt 1 title "(6,3) AB",\
     "minfoA.dat" using 1:5  with linespoints lc 2 pt 2 dt 3 title "(6,7) A",\
     "minfoA.dat" using 1:6  with linespoints lc 2 pt 2 dt 2 title "(6,7) B",\
     "minfoA.dat" using 1:7  with linespoints lc 2 pt 2 dt 1 title "(6,7) AB",\
     "minfoA.dat" using 1:8  with linespoints lc 3 pt 4 dt 3 title "(8,3) A",\
     "minfoA.dat" using 1:9  with linespoints lc 3 pt 4 dt 2 title "(8,3) B",\
     "minfoA.dat" using 1:10 with linespoints lc 3 pt 4 dt 1 title "(8,3) AB",\
     "minfoA.dat" using 1:11 with linespoints lc 4 pt 6 dt 3 title "(8,7) A",\
     "minfoA.dat" using 1:12 with linespoints lc 4 pt 6 dt 2 title "(8,7) B",\
     "minfoA.dat" using 1:13 with linespoints lc 4 pt 6 dt 1 title "(8,7) AB"