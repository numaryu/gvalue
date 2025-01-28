runname = "test_mixture"
outfile = "results_mix.eps"

nmedia = 2
ngeneration = 6

file_result = sprintf("%s.dat",runname)
file_degradation = sprintf("%s_degradation.dat",runname)
file_degradation_mix = sprintf("%s_degradation_mix.dat",runname)

set terminal postscript enhanced color eps "Times-Roman" 30 size 10in,10.5in
set output outfile


set multiplot layout 3, 2

set xlabel "Energy {/Times-Italic T} [eV]"
set ylabel "Stopping Power {/Times-Italic S}({/Times-Italic T}) [eV/cm]"

set log
set format x "10^{%L}"
set format y "10^{%L}"

set key top Left samplen 2 width 0

#set yrange [1.e-18:]

plot \
	file_result u 1:3 w l \
	lw 4 lc 1 \
	title "Mix", \
	for [i=1:nmedia] file_result index i-1 u 1:2 w l \
	lw 4 lc i+1 dt 2 \
	title sprintf("Medium %d",i)

set xlabel "Energy {/Times-Italic T} [eV]"
set ylabel "Slowing-down {/Times-Italic y}({/Times-Italic T}) [cm/eV]"

set log
set format x "10^{%L}"
set format y "10^{%L}"

set key bottom left Left samplen 2 width 0

#set yrange [1.e-18:]

plot \
	file_degradation_mix u 1:2 w lp \
	lw 4 lc 0 \
	title "Mix", \
	file_degradation u 1:2 w l \
	lw 4 lc 1 \
	title "Total", \
	for [i=1:ngeneration*nmedia] l = (i-1)/ngeneration file_degradation index l u 1:(column((i-1)%ngeneration+2)) w l \
	lw 4 lc (i-1)%ngeneration+2 dt l+1 \
	title sprintf("%d:%d",l+1,(i-1)%ngeneration+1, l+1)

set xlabel "Energy {/Times-Italic T} [eV]"
set ylabel "Total Cross Sectionn {/Times-Italic Q}({/Times-Italic T}) [cm^2]"

set log
set format x "10^{%L}"
set format y "10^{%L}"

set key bottom left Left samplen 2 width 0

#set yrange [1.e-18:]

plot \
	for [i=1:nmedia] file_result index i-1 u 1:5 w l \
	lw 4 lc 1 dt i \
	title sprintf("Ionization %d",i), \
	for [i=1:nmedia] file_result index i-1 u 1:6 w l \
	lw 4 lc 2 dt i\
	title sprintf("Singlet %d",i), \
	for [i=1:nmedia] file_result index i-1 u 1:7 w l \
	lw 4 lc 3 dt i \
	title sprintf("Triplet %d",i), \

set xlabel "Energy {/Times-Italic T} [eV]"
set ylabel "{/Times-Italic T} {/Times-Italic y} {/Times-Italic Q} [cm^3]"

unset log
set log x
set format x "10^{%L}"
set format y "% h"

set key top right Left samplen 2 width 0

plot \
	for [i=1:nmedia] file_result index i-1 u 1:8 w l \
	lw 4 lc 1 dt i \
	title sprintf("Ionization %d",i), \
	for [i=1:nmedia] file_result index i-1 u 1:9 w l \
	lw 4 lc 2 dt i \
	title sprintf("Singlet %d",i), \
	for [i=1:nmedia] file_result index i-1 u 1:10 w l \
	lw 4 lc 3 dt i \
	title sprintf("Triplet %d",i)

set xlabel "Energy {/Times-Italic T} [eV]"
set ylabel "Mean Free Path {/Symbol l}({/Times-Italic T}) [nm]"

set log x
set format x "10^{%L}"
set format y "% h"
set yrange [0.1:10]

set key top right Left samplen 2 width 0

factor_cm_to_nm=1.e7
plot \
	file_result u 1:($11*factor_cm_to_nm) w l \
	lw 4 lc 1 \
	title "Mix", \

set auto

set xlabel "Energy {/Times-Italic T} [eV]"
set ylabel "Range [nm]"

set log x
set format x "10^{%L}"
set format y "% h"
#set yrange [0.1:10]

set key top right Left samplen 2 width 0

factor_cm_to_nm=1.e7
plot \
	for [i=1:nmedia] file_result index i-1 u 1:($12*factor_cm_to_nm) w l \
	lw 4 lc 1 dt i \
	title sprintf("Ionization %d",i), \
	for [i=1:nmedia] file_result index i-1 u 1:($13*factor_cm_to_nm) w l \
	lw 4 lc 2 dt i \
	title sprintf("Singlet %d",i), \
	for [i=1:nmedia] file_result index i-1 u 1:($14*factor_cm_to_nm) w l \
	lw 4 lc 3 dt i \
	title sprintf("Triplet %d",i), \
