reset
set grid
set xlabel "P_d"
set ylabel "mutual info (per channel bit)"
set format x "%.3f"
set format y "%.2f"
set yrange [0.44:0.66]
set key left bottom
plot "rateC.dat" using 1:2 with linespoints lw 3 dt 1 title "code 1",\
     "rateC.dat" using 1:3 with linespoints lw 3 dt 2 title "code 2",\
     "rateC.dat" using 1:4 with linespoints lw 3 dt 3 title "code 3",\
     "rateC.dat" using 1:5 with linespoints lw 3 dt 4 title "code 4"