set xrange[0:5400]

dist = 8*1000*206264*150e6*1e5
factor = 4*3.1415*dist*dist
lumTomJy = 1e26/1.4e14

set xlabel 't [sec]'
set ylabel 'Flux at K band [mJy]'
plot 'lightCurveBlob.txt' u (5400-$1):($2/factor * lumTomJy) w l lw 3
