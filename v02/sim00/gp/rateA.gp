reset
set grid
set xlabel "P_i=P_d=P_s"
set ylabel "mutual information ({/Symbol n}=6)"
set format x "%.3f"
set format y "%.1f"
set yrange [2.7:3.9]
plot "rateA.dat" using 1:2 with linespoints lw 3 dt 1 title "code 1",\
     "rateA.dat" using 1:3 with linespoints lw 3 dt 2 title "code 2",\
     "rateA.dat" using 1:4 with linespoints lw 3 dt 3 title "code 3"