reset
set grid
set xlabel "P_{id}"
set ylabel "R"
set yrange [0.4:0.8]
set format x "%.3f"
set format y "%.2f"
plot "../data/MK3_1_0_010.txt"  using 1:4 with linespoints title "MK3.1.0",\
     "../data/MK4_1_0_010.txt"  using 1:4 with linespoints title "MK4.1.0",\
     "../data/MK4_2_2_010.txt"  using 1:4 with linespoints title "MK4.2.2",\
     "../data/MK5_1_0_010.txt"  using 1:4 with linespoints title "MK5.1.0",\
     "../data/MK6_1_0_010.txt"  using 1:4 with linespoints title "MK6.1.0",\
     "../data/MK6_2_2_010.txt"  using 1:4 with linespoints title "MK6.2.2",\
     "../data/MK7_1_0_010.txt"  using 1:4 with linespoints title "MK7.1.0",\
     "../data/MK8_1_0_010.txt"  using 1:4 with linespoints title "MK8.1.0",\
     "../data/MK8_2_2_010.txt"  using 1:4 with linespoints title "MK8.2.2",\
     "../data/MK10_2_2_010.txt" using 1:4 with linespoints title "MK10.2.2"