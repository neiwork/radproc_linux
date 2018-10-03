	set xrange [7.5:26]
	set yrange [32:48]
	set xlabel "frecuency [Hz]"
	set ylabel "vLv [erg s^-1]

	plot 'lum.txt' u (log10($1)):(log10($2)) t 'Sync' w l lw 2, 'lum.txt' u (log10($1)):(log10($3)) t 'Bremss' w l lw 2, 'lum.txt' u (log10($1)):(log10($4)) t 'IC_{in}' w l lw 2, 'lum.txt' u (log10($1)):(log10($5)) t 'IC_{Out}' w l lw 2, 'lum.txt' u (log10($1)):(log10($6)) t 'Tot' w l lw 2
