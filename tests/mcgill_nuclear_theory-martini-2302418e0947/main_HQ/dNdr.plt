set terminal postscript eps enh color
set output "dNdr.eps"
set xlabel "r [fm]"
set logscale y
plot "dN4pir2dr_ccbar_beforeEvolution_RHIC_AuAu_0t20_evolved_1.dat" with lines lw 4 t "Before evolution", \
"dN4pir2dr_ccbar_RHIC_AuAu_0t20_evolved_1.dat" with lines lw 4 t "After evolution"
q