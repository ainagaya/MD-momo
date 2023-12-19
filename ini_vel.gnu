# BASE PLOTTER
# A partir d'un script creat al Novembre del 2019 (Autora: Aina Gaya Àvila)

# Indicacions inicials, tamany de la figura, tipus de lletra i nom del fitxer de sortida
set terminal pngcairo size 820,500 enhanced font 'Verdana,13'
set output "ini_vel_x.png"
set style data histogram
set style fill solid

#-------------------------------------
#------Definició d'estils (més dels necessaris)
set style line 1 lc rgb '#8b1a0e' pt 1 ps 1 lt 1 lw 1 dt 2 # --- red
set style line 2 lc rgb '#5e9c36' pt 2 ps 1 lt 2 lw 1 # --- green
set style line 3 lc rgb '#1c6ee8' pt 3 ps 1 lt 3 lw 1 # --- blue
set style line 4 lc rgb '#4F1271' pt 3 ps 1 lt 3 lw 1 # --- indigo
set style line 5 lc rgb '#1D1E18' pt 3 ps 1 lt 3 lw 1 # --- black


#--------------------------------------
#-- X-Axis
set xlabel "Number of occurrences"      # Etiqueta de l'eix X

#--------------------------------------
#-- Y-Axis
set ylabel "Velocity"

#--------------------------------------
#-- Labels and Title
set title 'Initial Velocities'

#--------------------------------
#--- PLOT

file = "vel_ini.dat"

n=100 #number of intervals
max=15. #max value
min=-15. #min value
width=(max-min)/n #interval width
#function used to map a value to the intervals
hist(x,width)=width*floor(x/width)+width/2.0
set boxwidth width*0.9
set style fill solid 0.5 # fill style

#count and plot
plot file u (hist($1,width)):(1.0) smooth freq w boxes ls 1 t 'x'

set output "ini_vel_y.png"
plot file u (hist($2,width)):(1.0) smooth freq w boxes ls 2 t 'y'

set output "ini_vel_z.png"
plot file u (hist($3,width)):(1.0) smooth freq w boxes ls 3 t 'z'

n=30 #number of intervals
max=30. #max value
min=-30. #min value
width=(max-min)/n #interval width
#function used to map a value to the intervals
hist(x,width)=width*floor(x/width)+width/2.0
set boxwidth width*0.9
set style fill solid 0.5 # fill style

# ---------------------

T = 85.33

# Define the Gaussian distribution
f(x) = 300*(1./(2.*3.1416*T))**(1./2.)*exp(-(x**2)/(2.*T))
g(x) = 300*(1./(2.*3.1416*T))**(3./2.)*4.*3.1416*(x+2)**2*exp(-(x+2)**2/(2.*T))

# Set the output file names
set output "fin_vel_x_verlet.png"
file = "vel_fin_Verlet.dat"

set xlabel "v_x"
set ylabel "occurrences"

# Plot the normalized histogram and the Gaussian distribution for x
plot file i 1 u (hist($1,width)) : (1.0) smooth freq w boxes ls 1 t 'x', f(x) t'M-B' ls 5

set output "fin_vel_y_verlet.png"
set xlabel "v_y"
# Plot the normalized histogram and the Gaussian distribution for y
plot file i 1 u (hist($2,width)) : (1.0) smooth freq w boxes ls 2 t 'y', f(x) t'M-B' ls 5

set output "fin_vel_z_verlet.png"
set xlabel "v_z"
# Plot the normalized histogram and the Gaussian distribution for z
plot file i 1 u (hist($3,width)) : (1.0) smooth freq w boxes ls 3 t 'z', f(x) t'M-B' ls 5

set output "fin_vel_module_verlet.png"
set xlabel "v"
# Plot the normalized histogram and the Gaussian distribution for the module
plot file i 1 u (hist($4,width)) : (1.0) smooth freq w boxes ls 4 t 'v', g(x) t'M-B' ls 5

