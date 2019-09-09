set xrange[1000:6500]

dist = 8000*206264*150e11
factor = 4*3.1415*dist*dist
lumTomJy = 1.0e26/1.4e14

set xlabel 't [sec]'
set ylabel 'Flux at K band [mJy]'
plot 'lightCurveBlob.txt' u 1:($2/factor * lumTomJy) w l lw 3
