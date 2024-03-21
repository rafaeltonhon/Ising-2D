set xl ' $b$'
set yl ' $\sigma_b^2$'
set key box lw 2 lc black
set grid lw 3 lc 'gray' lt 1
set ytics format '%.4f'
set key spacing 1.2
set key Left
set parametric
set yr[:0.0006]
set xr[0:80]
set tr[0:100]
set terminal epslatex color colortex standalone size 7,6 font ',20'
set output 'bining.tex'
f(t)=a+b/t
g(t)=c+d/t
h(t)=e+f1/t
i(t)=y+z/t

#fit [10:100] f(t) 'bining.dat' via a, b
fit [20:100] g(t) 'bining.dat' via c, d
fit [50:100] h(t) 'bining.dat' via e, f1
fit [0:100] i(t) 'bining.dat' via y, z

plot 'bining.dat' lw 3 pt 4 ps 1.5 lc 'black' title ' $\sigma_b^2$',\
    t,i(t) lw 8 lc 'red' t ' $b_0=0$',\
    t,g(t) lw 8 lc 'blue' t ' $b_0=20$',\
    t,h(t) lw 8 lc 12 t ' $b_0=50$'
set output