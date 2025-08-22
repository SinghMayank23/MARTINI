set terminal postscript eps enh color
set output "dNdy.eps"
set xlabel "y"
plot "dNdy_ccbar_beforeEvolution_RHIC_AuAu_0t20_evolved_1.dat" with lines lw 4 t "Before evolution", \
"dNdy_ccbar_RHIC_AuAu_0t20_evolved_1.dat" with lines lw 4 t "After evolution"
q