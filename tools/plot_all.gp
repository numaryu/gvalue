runname = "test"
outfile = "results.eps"

file_result = sprintf("%s.dat",runname)
file_degradation = sprintf("%s_degradation.dat",runname)

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
	file_result u 1:2 w l \
	lw 4 lc 1 \
	title "", \

set xlabel "Energy {/Times-Italic T} [eV]"
set ylabel "Slowing-down {/Times-Italic y}({/Times-Italic T}) [cm/eV]"

set log
set format x "10^{%L}"
set format y "10^{%L}"

set key bottom left Left samplen 2 width 0

#set yrange [1.e-18:]

plot \
	file_degradation u 1:2 w l \
	lw 4 lc 1 \
	title "Total", \
	file_degradation u 1:3 w l \
	lw 4 lc 2 dt 2 \
	title "1", \
	file_degradation u 1:4 w l \
	lw 4 lc 3 dt 2 \
	title "2", \
	file_degradation u 1:5 w l \
	lw 4 lc 4 dt 2 \
	title "3", \
	file_degradation u 1:6 w l \
	lw 4 lc 5 dt 2 \
	title "4", \
	file_degradation u 1:7 w l \
	lw 4 lc 6 dt 2 \
	title "5", \
	file_degradation u 1:8 w l \
	lw 4 lc 7 dt 2\
	title "6"

set xlabel "Energy {/Times-Italic T} [eV]"
set ylabel "Total Cross Sectionn {/Times-Italic Q}({/Times-Italic T}) [cm^2]"

set log
set format x "10^{%L}"
set format y "10^{%L}"

set key bottom left Left samplen 2 width 0

#set yrange [1.e-18:]

plot \
	file_result u 1:4 w l \
	lw 4 lc 1 \
	title "Ionization", \
	file_result u 1:5 w l \
	lw 4 lc 2 \
	title "Singlet", \
	file_result u 1:6 w l \
	lw 4 lc 3 \
	title "Triplet", \

set xlabel "Energy {/Times-Italic T} [eV]"
set ylabel "{/Times-Italic T} {/Times-Italic y} {/Times-Italic Q} [cm^3]"

unset log
set log x
set format x "10^{%L}"
set format y "% h"

set key top right Left samplen 2 width 0

plot \
	file_result u 1:7 w l \
	lw 4 lc 1 \
	title "Ionization", \
	file_result u 1:8 w l \
	lw 4 lc 2 \
	title "Singlet", \
	file_result u 1:9 w l \
	lw 4 lc 3 \
	title "Triplet"

set xlabel "Energy {/Times-Italic T} [eV]"
set ylabel "Mean Free Path {/Symbol l}({/Times-Italic T}) [nm]"

set log x
set format x "10^{%L}"
set format y "% h"
set yrange [0.1:10]

set key top right Left samplen 2 width 0

factor_cm_to_nm=1.e7
plot \
	file_result u 1:($10*factor_cm_to_nm) w l \
	lw 4 lc 1 \
	title "", \

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
	file_result u 1:($11*factor_cm_to_nm) w l \
	lw 4 lc 1 \
	title "Ionization", \
	file_result u 1:($12*factor_cm_to_nm) w l \
	lw 4 lc 2 \
	title "Singlet", \
	file_result u 1:($13*factor_cm_to_nm) w l \
	lw 4 lc 3 \
	title "Triplet", \
