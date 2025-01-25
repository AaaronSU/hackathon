set terminal pngcairo enhanced


set xlabel 'Number of cores (log2)'
set ylabel 'Time (s)'
set y2label 'Mean results with Error bar'
set logscale x 2
set y2range [5: 6]
set y2tics

set style data linespoints

set output 'ss7.png'
set title 'Strong Scaling c7g 100000 1000000'

plot 'ss7.txt' using 1:2 with linespoints title 'Time (s)', \
     'ss7.txt' using 1:3 with linespoints axes x1y2 title 'Mean results', \
     'ss7.txt' using 1:3:4 with errorbars axes x1y2 title 'Error bar', \
     3288.152810 / x title 'strong scaling reference' with lines lc 'red' linewidth 2 

set output 'ss8.png'
set title 'Strong Scaling c8g 100000 1000000'

plot 'ss8.txt' using 1:2 with linespoints title 'Time (s)', \
     'ss8.txt' using 1:3 with linespoints axes x1y2 title 'Mean results', \
     'ss8.txt' using 1:3:4 with errorbars axes x1y2 title 'Error bar', \
     2891.242000 / x  title 'strong scaling reference' with lines lc 'red' linewidth 2

set output 'ws7.png'
set title 'Weak Scaling c7g 10000*Cores 1000000'
set yrange [200: 400]

plot 'ws7.txt' using 1:2 with linespoints title 'Time (s)', \
     'ws7.txt' using 1:3 with linespoints axes x1y2 title 'Mean results', \
     'ws7.txt' using 1:3:4 with errorbars axes x1y2 title 'Error bar', \
     331.844850 title 'weak scaling reference' with lines lc 'red' linewidth 2

set output 'ws8.png'
set title 'Weak Scaling c8g 10000*Cores 1000000'
set yrange [200: 400]

plot 'ws8.txt' using 1:2 with linespoints title 'Time (s)', \
     'ws8.txt' using 1:3 with linespoints axes x1y2 title 'Mean results', \
     'ws8.txt' using 1:3:4 with errorbars axes x1y2 title 'Error bar', \
     290.039931 title 'weak scaling reference' with lines lc 'red' linewidth 2