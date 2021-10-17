reset
#set terminal pdf
set term postscript color eps font "Times-Roman,20"
set tics font "Times-Roman,20"
set output 'scal_ssna303.ps'
set xlabel 'Nb of cores'
set ylabel 'Time (s)'
set format y "10^{%L}" 
set xtics 2 
f(x) = x**(-1)*4000

set logscale


set title "\
Notched beam case, 1.5M unknowns \n\
Strong scaling, Execution time (8 time steps)\n\
on mesocentre Aix-Marseille Univ." offset -1

plot [16:256] 'petsc.txt' u 1:2 w lp t "PETSc with rc\\_ex10p setting" ls 7  lw 2 lc rgb 'blue',\
'mumps.txt' u 1:2 w lp t "MUMPS direct solver" ls 7  lw 2 lc rgb 'green',\
f(x) t "ideal speedup" ls 3 lc rgb 'black' dt 2 lw 3,\
'hypre.txt' u 1:2 w lp t "HypreFGMRES solver + BoomerAMG precond" ls 7  lw 2 lc rgb 'red'

# 'sca u 1:4 w lp t "Lecture maillage" ls 7 ps 2 lw 2 lc rgb 'purple', \
#'sca u 1:6 w lp t "Assemblage" ls 7 ps 2 lw 2 lc rgb 'yellow',\
#'sca u 1:7 w lp t "Resolution"  ls 7 ps 2 lw 2 lc rgb 'green'
