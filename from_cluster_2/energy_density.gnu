#BASE PLOTTER
#A partir d'un script creat al Novembre del 2019 (Autora: Aina Gaya Àvila)


#Indicacions inicials, tamany de la figura, tipus de lletra i nom del fitxer de sortida
set terminal pngcairo size 820,500 enhanced font 'Verdana,13'
set output "pot_energy_density.png"

set zero 1.e-20               # Agafa com zero quantitats inferiors
#set encoding  iso_8859_15     # Permet l'ús d'accents als  MAS 
#set encoding  iso_8859_1      # Permet l'ús d'accents als PC's    
#set encoding  CP1250          # Permet l'ús d'accents als MS-Windows
#set locale "es_ES.utf8"

#-------------------------------------------------------------------
#------ Marges, graella 

#set lmargin at screen 0.15 
#set rmargin at screen 0.9
#set bmargin at screen 0.15
#set tmargin at screen 0.8 

#Treu borde de dalt i dret. La resta els fa grisos.
#set style line 11 lc rgb '#808080' lt 1
#set border 3 back ls 11
#set tics nomirror
#Defineix la graella, també grisa
#set style line 12 lc rgb '#808080' lt 0 lw 1
#set grid back ls 12
 
#-------------------------------------
#------ Llegenda 

set key center right                    # Llegenda al costat esquerre 
#set key font "Arial, 20"        # Font de la legenda

#set key spacing 1.5 R font  "Arial, 18" # Espaiat de línies i tamany font a la legenda

#set key outside            # Llegenda en aquestes coordenades
#set key box lt -1 lw 1           # Recepte per posar una caixa a la legenda
#set key spacing 1.5 R            # Espaiat de línies

#-------------------------------------
#------Definició d'estils (més dels necessaris)

set style line 1 lc rgb '#8b1a0e' pt 1 ps 1 lt 1 lw 1 dt 2 # --- red
set style line 2 lc rgb '#5e9c36' pt 2 ps 1 lt 2 lw 1 # --- green
set style line 3 lc rgb '#1c6ee8' pt 3 ps 1 lt 3 lw 1 # --- blue
set style line 4 lc rgb '#ffa500' pt 3 ps 1 lt 4 lw 1 # --- orange
set style line 5 lc rgb '#9370DB' pt 5 ps 1 lt 5 lw 1 # --- purple
set style line 6 lc rgb '#FA8072' pt 6 ps 1 lt 6 lw 1 # --- salmon
set style line 7 lc rgb '#738EFA' pt 7 ps 1 lt 7 lw 1 # --- light blue
set style line 8 lc rgb '#FA7373' pt 1 ps 1 lt 8 lw 1 # --- light red
set style line 9 lc rgb '#73FA7C' pt 9 ps 1 lt 9 lw 1 # --- light green


#--------------------------------------
#-- X-Axis

set xlabel "Density (g/cm^3)"      # Etiqueta de l'eix X
#set xrange [0.:220.]          # Limits de l'eix X

#set xtics 10e-11          
#set mxtics 5                 # Numero de subtics per cada xtic
#set xtics border mirror 10

#set grid xtics                 # Linies de grid X fix
#set xtics  font "Arial, 20"    # Font per les etiquetes situades als xtics
#set xlabel font "Arial, 20"     
#set logscale x                 # Activa l'escala logaritmica

#--------------------------------------
#-- Y-Axis

#set ytics  font "Arial, 20"
#set ylabel font "Arial, 20"
#set yrange [0:100]
#set ytics 5.e5,1.e5,1.e6
#set mytics 5
#set ytics border mirror 1

set ylabel "E (kJ/mol)"
#set grid ytics
#set format y "10^%L"        # Les etiquetes logaritmiques són millors.
#set logscale y

#--------------------------------------
#-- Labels and Title

#set label 1 "bla bla " at screen 0.5,0.69 center font "Helvetica, 24"

set title "Energy vs. density in a Lennard-Jones liquid at k_B T=1.2"
#set title font "Arial, 16"

#set label "S^{3}_{N} = Número de porcions amb N-2 pizzes" at screen 0.35, 0.83 font "Arial,18"

#--------------------------------
#--- PLOT

data_file = "results.dat"

epsilon = 0.998
sigma = 3.4e-8
m = 40
Avog = 6.023e23
part = 125

to_pascal = 10**33 * part/Avog
to_pascal = part/Avog

plot data_file using ($1*m/sigma**3/Avog*part):($3*epsilon):($4*epsilon) with yerrorbars title 'Kinetic', \
     data_file using ($1*m/sigma**3/Avog*part):($5*epsilon):($6*epsilon) with yerrorbars title 'Potential', \
     data_file using ($1*m/sigma**3/Avog*part):($7*epsilon):($8*epsilon) with yerrorbars title 'Total'
     
     
set output "preassure.png"
set title "Preassure vs. density in a Lennard-Jones liquid at k_B T=1.2"
set ylabel "P (Pa)"

set xzeroaxis

plot data_file using ($1*m/sigma**3/Avog*part):(($9*epsilon/sigma**3)*to_pascal):(($10*epsilon/sigma**3)*to_pascal) with yerrorbars t''

     
     
     
     
     
