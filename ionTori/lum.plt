	set xrange [5:25]
	set yrange [25:45]
	set xlabel "frecuency [Hz]"
	set ylabel "vLv [erg s^-1]

	plot 'lum.txt' u (log10($1)):(log10($2)) t 'Sync' w l lw 2, 'lum.txt' u (log10($1)):(log10($3)) t 'Bremss' w l lw 2, 'lum.txt' u (log10($1)):(log10($6)) t 'IC_B' w l lw 2, 'lum.txt' u (log10($1)):(log10($4)) t 'IC_{inS}' w l lw 2, 'lum.txt' u (log10($1)):(log10($5)) t 'IC_{inB}' w l lw 2,
