set terminal wxt
set grid

set tics out
set y2tics

x0=NaN
y0=NaN

plot '../intercept_distance.txt' using 1:2 with lines title 'd' , \
    '' using 1:3 with lines title 'dd' axes x1y2, \
    '' using 1:4 with lines title 'ddd' axes x1y2,\
    '' using 1:5 with lines title 'tgt',\
    '../intercept_steps.txt' using 2:3:1 title 'steps' with labels

#'' using 1:6 with lines title 'newton' axes x1y2,\
#'' using 1:7 with lines title 'halley' axes x1y2,\
#'' using 1:8 with lines title 'laguerre' axes x1y2,\
#'' using (dx=$1-x0,x0=$1,$1):(dy=$2-y0,y0=$2,dy/dx) with lines title 'derivative v',\
#'' using (dx=$1-x0,x0=$1,$1):(dy=$3-y0,y0=$3,dy/dx) with lines title 'derivative a'


#'' using 1:(-$5) with lines title 'acceleration', \
