reset 

set xtics offset 0,-1,0
unset key

set view equal xy 
set view 0,0,2

set terminal png size 1024,768 crop enhanced

do for [i=0:5]{

output_file = sprintf('T_%d_animation.png',i)

set output output_file

splot 'T'.i.'.dat' u 1:2:6 notitle w pm3d 
}

set output

reset

set terminal pop

