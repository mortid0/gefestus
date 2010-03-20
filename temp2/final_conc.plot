set size 1,1
set terminal postscript enhanced
set title 'Counterflow mass'
set yrange[0:1.1]
set xrange[100:7200]
#set format x ''
#set format y ''
set border 3
set label 'Diam, m' at 6000,0.1
set label ' ~M{.6\~}' at 200,1.1
set xtics nomirror
set ytics nomirror
unset mx2tics 
unset my2tics
unset x2tics 
unset y2tics
#unset xtics
#unset ytics
#unset xlabel 
#unset ylabel
set key 2500,0.8
set output 'final_conc.eps'
plot 'plot_num_0.dat' u 1:(($9-81.372)/1016.628) with l lw 2 ti 'N_0 = 0.63 E-15, m^{-3}', 'plot3.dat' u 1:(($9-98.093)/1197.907)  every 2:2 w l lw 2 ti 'N_0 = 1.26 E-15, m^{-3}', 'plot2.dat' u 1:(($9-41.443)/435.607) w l lw 2 ti 'N_0 = 0.32 E-15, m^{-3}'