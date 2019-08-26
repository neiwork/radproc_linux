set xrange [8:21]
set yrange [-15:-10]
dist = 8*1000*206265*150e6*1.0e5
factor = 4.0*3.1415*dist*dist

set xlabel 'Freq [Hz]'
set ylabel 'nuFnu [erg cm^{-2} s^{-1}]'

plot 'lum0.txt' u (log10($1)):(log10(($3+$4+$5)/factor)) w l lw 2, 'lum1.txt' u (log10($1)):(log10(($3+$4+$5)/factor)) w l lw 2, 'lum2.txt' u (log10($1)):(log10(($3+$4+$5)/factor)) w l lw 2, 'lum3.txt' u (log10($1)):(log10(($3+$4+$5)/factor)) w l lw 2, 'lum4.txt' u (log10($1)):(log10(($3+$4+$5)/factor)) w l lw 2, 'lum5.txt' u (log10($1)):(log10(($3+$4+$5)/factor)) w l lw 2

