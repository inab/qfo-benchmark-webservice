set terminal png medium size 800,800 font "FreeSerif,12" enhanced
set output '<DVAR name="outfile">'
set format y "%4.3f"
set xlabel "<DVAR name="xlabel">"
set ylabel "<DVAR name="ylabel">"
set key outside bottom center
set style line 200 lc rgb "#333333" lw 2 pt 1
set style arrow 1 head filled size screen 0.015,25,55 ls 200
set arrow from graph 0.90,0.85 to graph 0.95,0.95 as 1
set label "better" at graph 0.92,0.92 right


plot \
<DLOOP name="plot"> \
  '<DVAR name="datafile">' i <DLVAR name="index"> u 2:3:4:5 tit '<DLVAR name="title">' w xyerror ls <DLVAR name="style">,\
</DLOOP> \
1/0 notit
