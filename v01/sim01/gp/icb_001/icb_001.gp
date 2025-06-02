reset
set grid
set xlabel "P_{id}"
set ylabel "R"
set yrange [0.4:0.8]
set format x "%.3f"
set format y "%.2f"
plot "../data/CB3_4_2_001.txt"  using 1:4 with linespoints title "CB3.4.2",\
     "../data/CB3_5_1_001.txt"  using 1:4 with linespoints title "CB3.5.1",\
     "../data/CB4_5_2_001.txt"  using 1:4 with linespoints title "CB4.5.2",\
     "../data/CB4_6_1a_001.txt" using 1:4 with linespoints title "CB4.6.1a",\
     "../data/CB4_6_1b_001.txt" using 1:4 with linespoints title "CB4.6.1b",\
     "../data/CB4_7_1_001.txt"  using 1:4 with linespoints title "CB4.7.1b",\
     "../data/CB5_7_2_001.txt"  using 1:4 with linespoints title "CB5.7.2",\
     "../data/CB5_a_2_001.txt"  using 1:4 with linespoints title "CB5.a.2",\
     "../data/CB5_8_2_001.txt"  using 1:4 with linespoints title "CB5.8.2",\
     "../data/CB5_9_1_001.txt"  using 1:4 with linespoints title "CB5.9.1"