#Reset gnuplot:
reset session

#Titel plot:
set title "Hele mooie output" 

#Set view { <rot_x>{,{<rot_z>}{,{<scale>}{,<scale_z>}}} | map }
set view 180, 0, 2, 1

#Hide grid lines:
set hidden3d

#Lines instead of grid points:
set style data lines

#Contour lines:
set contour
 # Position legend countours
 set key outside right bottom        # {left | right | center} {top | bottom | center}


#Range of axis:
set xrange [ 0.0 : 2.0 ] noreverse nowriteback
set yrange [ 0 : 0.2 ] noreverse nowriteback
#set zrange [ 0 : 0.8 ] noreverse nowriteback

#Set equal scale along all three axis
set view equal xyz

#Labels:
set xlabel "X axis" 
set xlabel  offset character -3, -2, 0 font "" textcolor lt -1 norotate
set ylabel "Y axis" 
set ylabel  offset character 3, -2, 0 font "" textcolor lt -1 norotate
set zlabel "Z axis" 
set zlabel  offset character -3, 0, 0 font "" textcolor lt -1 norotate

#Plot comluns form .dat file:

#splot 'output.dat' using 1:2:12 with pm3d

scale = 0.15
#plot 'output.dat' using 1:2:($3*scale):($4*scale) with vectors head filled lw 2


#Check sneltoetsen van gnuplot:
#show bind

### Start multiplot 
set multiplot layout 1,1 rowsfirst

#set label 1 'ugrid' at graph 0.0,0.0,0.0 font ',16'


splot 'output.dat' using 1:2:3 with pm3d # --- GRAPH: ugrid

#splot 'output.dat' using 1:2:4 with pm3d # --- GRAPH: vgrid

#splot 'output.dat' using 1:2:5 with pm3d # --- GRAPH: p

#splot 'output.dat' using 1:2:6 with pm3d # --- GRAPH: T

#splot 'output.dat' using 1:2:7 with pm3d # --- GRAPH: rho

#splot 'output.dat' using 1:2:8 with pm3d # --- GRAPH: mu

#splot 'output.dat' using 1:2:9 with pm3d # --- GRAPH: Gamma

#splot 'output.dat' using 1:2:10 with pm3d # --- GRAPH: k

#splot 'output.dat' using 1:2:11 with pm3d # --- GRAPH: eps

#splot 'output.dat' using 1:2:12 with pm3d # --- GRAPH: Tplus_u

#splot 'output.dat' using 1:2:13 with pm3d # --- GRAPH: Tplus_v

#splot 'output.dat' using 1:2:14 with pm3d # --- GRAPH: yplus_u

#splot 'output.dat' using 1:2:15 with pm3d # --- GRAPH: yplus_v

#splot 'output.dat' using 1:2:16 with pm3d # --- GRAPH: uplus_u

#splot 'output.dat' using 1:2:17 with pm3d # --- GRAPH: uplus_v

#splot 'vort.dat' using 1:2:3 with pm3d # --- GRAPH: x[I], y[J], vorticity

#splot 'str.dat' using 1:2:3 with pm3d # --- GRAPH: x[I], y[J], stream

#splot 'velu.dat' using 1:2:3 with pm3d # --- GRAPH: x_u[i], y[J], u[i][J]

#splot 'velv.dat' using 1:2:3 with pm3d # --- GRAPH: x[I], y_v[j], v[I][j]

unset multiplot
### End multiplot

## Create animation:
#set style data lines
#do for [i=1:20] { splot sprintf('data%d.dat', i) using 4:5:6; pause 0.5 }

