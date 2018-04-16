reset 

set xtics offset 0,-1,0
unset key

set view equal xy 
set view 0,0,2

set cbrange [250:320] # Temp
#set cbrange [-1:4]   # Vel_u

set terminal png size 1024,768 crop enhanced

do for [i=0:140]{

output_file = sprintf('T_%d_animation.png',i)

set output output_file

splot 'T'.i.'.dat' u 1:2:6 notitle w pm3d 
}

set output

reset

set terminal pop

