f(x) = a+b*x+c*x*x+d*x*x*x+e*x*x*x*x+f*x*x*x*x*x+g*x*x*x*x*x*x
fit f(x) 'newTempElectrons.txt' u 1:2 via a,b,c,d,e,f,g

set table 'output.txt'
set samples 100
set xrange [0.01:2.01]

plot f(x)
unset table

set yrange [8:10]
plot 'fields.dat' w l lw 5, 'newTempElectrons.txt', f(x), 'fields_original.dat' w l lw 3
