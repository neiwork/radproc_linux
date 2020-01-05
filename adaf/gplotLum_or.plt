set xrange [12:22]
set yrange [35:45]

plot 'lum.txt' u (log10($1)):(log10($3)) w l lw 3, 'lum.txt' u (log10($1)):(log10($4)) w l lw 4, 'lum.txt' u (log10($1)):(log10($5)) w l lw 4, 'lum_or.txt' u (log10($1)):(log10($3)) w l lw 3, 'lum_or.txt' u (log10($1)):(log10($4)) w l lw 4, 'lum_or.txt' u (log10($1)):(log10($5)) w l lw 4, 'lum.txt' u (log10($1)):(log10($7)) w l lw 2, 'lum.txt' u (log10($1)):(log10($8)) w l lw 2, 'lum.txt' u (log10($1)):(log10($6)) w  l lw 2, 'lum_or.txt' u (log10($1)):(log10($6)) w l lw 2

