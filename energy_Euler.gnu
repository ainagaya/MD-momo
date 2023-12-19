# BASE PLOTTER
# A partir d'un script creat al Novembre del 2019 (Autora: Aina Gaya Àvila)

# Indicacions inicials, tamany de la figura, tipus de lletra i nom del fitxer de sortida
set terminal pngcairo size 820,500 enhanced font 'Verdana,13'
set output "energy_euler_1.png"

set zero 1.e-20               # Agafa com zero quantitats inferiors

#-------------------------------------
#------ Llegenda 
set key center right  

#-------------------------------------
#------Definició d'estils (més dels necessaris)
set style line 1 lc rgb '#702632' pt 1 ps 1 lt 1 lw 3 # --- red
set style line 2 lc rgb '#28AFB0' pt 2 ps 1 lt 2 lw 3 # --- rose
set style line 3 lc rgb '#AF9AB2' pt 3 ps 1 lt 3 lw 3 # --- blue
set style line 4 lc rgb '#F46197' pt 3 ps 1 lt 3 lw 3 # --- rosa chicle


#--------------------------------------
#-- X-Axis
set xlabel "t'"      # Etiqueta de l'eix X

#--------------------------------------
#-- Y-Axis
set ylabel "E'"

#--------------------------------------
#-- Labels and Title
set title "\Delta t' = 10^{-3}"

#--------------------------------
#--- PLOT
file = "energy_euler.dat"

plot file i 0 using 1:2 with lines linestyle 1 t 'Potential energy', \
     file i 0 using 1:3 with lines linestyle 2 t 'Kinetic energy', \
     file i 0 using 1:4 with lines linestyle 3 t 'Total energy'

set output "energy_euler_2.png"
set title "\Delta t' = 10^{-4}"

plot file i 1 using 1:2 with lines linestyle 1 t 'Potential energy', \
     file i 1 using 1:3 with lines linestyle 2 t 'Kinetic energy', \
     file i 1 using 1:4 with lines linestyle 3 t 'Total energy'

set output "energy_euler_3.png"
set title "\Delta t' = 10^{-5}"

plot file i 2 using 1:2 with lines linestyle 1 t 'Potential energy', \
     file i 2 using 1:3 with lines linestyle 2 t 'Kinetic energy', \
     file i 2 using 1:4 with lines linestyle 3 t 'Total energy'

# MOMENTUM PLOTS

unset key
set format y "%15.13f"

set output "momentum_euler_1.png"
set title "\Delta t' = 10^{-3}"
set ylabel "p'"

plot file i 0 using 1:5 with lines linestyle 4 t 'Momentum'

set output "momentum_euler_2.png"
set title "\Delta t' = 10^{-4}"

plot file i 1 using 1:5 with lines linestyle 4 t 'Momentum'

set output "momentum_euler_3.png"
set title "\Delta t' = 10^{-5}"

plot file i 2 using 1:5 with lines linestyle 4 t 'Momentum'

