set xrange [12:22]
set yrange [35:45]

plot 'lumThermal.dat' u (log10($1)):(log10($3)) w l lw 3, 'lumThermal.dat' u (log10($1)):(log10($4)) w l lw 3, 'lumThermal.dat' u (log10($1)):(log10($5)) w l lw 3, 'lumThermal.dat' u (log10($1)):(log10($9)) w l lw 4
