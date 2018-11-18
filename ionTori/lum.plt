	set xrange [7.5:26]
	set yrange [32:48]
	set xlabel "energy [eV]"
	set ylabel "vLv [erg s^-1]

	plot 'lum.txt' u (log10($2)):(log10($3)) t 'Sync' w l lw 2, 'lum.txt' u (log10($2)):(log10($4)) t 'Bremss' w l lw 2, 'lum.txt' u (log10($2)):(log10($6)) t 'IC_{Out}' w l lw 2, 'lum.txt' u (log10($2)):(log10($7)) t 'pp' w l lw 2, 'lum.txt' u (log10($2)):(log10($8)) t 'Tot' w l lw 2, 'lum.txt' u (log10($2)):(log10($5)) t 'IC_in' w l lw 2
