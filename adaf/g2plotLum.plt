set xrange [12:22]
set yrange [18:25]

plot 'lum.txt' u (log10($1)):(log10($3/$1)) w l lw 3, 'lum.txt' u (log10($1)):(log10($4/$1)) w l lw 4, 'lum.txt' u (log10($1)):(log10($5/$1)) w l lw 4
