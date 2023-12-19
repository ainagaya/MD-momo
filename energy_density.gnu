# BASE PLOTTER
# A partir d'un script creat al Novembre del 2019 (Autora: Aina Gaya Àvila)

# Indicacions inicials, tamany de la figura, tipus de lletra i nom del fitxer de sortida
set terminal pngcairo size 820,500 enhanced font 'Verdana,13'
set output "pot_energy_density.png"

set zero 1.e-20               # Agafa com zero quantitats inferiors

#-------------------------------------
#------ Llegenda 
set key center right                    # Llegenda al costat esquerre 

#-------------------------------------
#------Definició d'estils (més dels necessaris)
set style line 1 lc rgb '#702632' pt 1 ps 1 lt 1 lw 3 # --- red
set style line 2 lc rgb '#28AFB0' pt 2 ps 1 lt 2 lw 3 # --- rose
set style line 3 lc rgb '#AF9AB2' pt 3 ps 1 lt 3 lw 3 # --- blue
set style line 4 lc rgb '#F46197' pt 3 ps 1 lt 3 lw 3 # --- rosa chicle


#--------------------------------------
#-- X-Axis
set xlabel "Density (g/cm^3)"      # Etiqueta de l'eix X

#--------------------------------------
#-- Y-Axis
set ylabel "E (kJ/mol)"

#--------------------------------------
#-- Labels and Title
set title "Energy vs. density in a Lennard-Jones liquid at k_B T=1.2"

#--------------------------------
#--- PLOT
data_file = "results.dat"

epsilon = 0.998
sigma = 3.4e-10
m = 40
Avog = 6.023e23
part = 125

to_pascal = 10e3 * part/Avog

plot data_file using ($1*m/(sigma*100)**3/Avog*part):($3*epsilon):($4*epsilon) with yerrorbars linestyle 1 title 'Kinetic', \
     data_file using ($1*m/(sigma*100)**3/Avog*part):($5*epsilon):($6*epsilon) with yerrorbars linestyle 2 title 'Potential', \
     data_file using ($1*m/(sigma*100)**3/Avog*part):($7*epsilon):($8*epsilon) with yerrorbars linestyle 3 title 'Total'

set output "preassure.png"
set title "Pressure vs. density in a Lennard-Jones liquid at k_B T=1.2"
set ylabel "P (Pa)"

set xzeroaxis

plot data_file using ($1*m/(sigma*100)**3/Avog*part):(($9*epsilon/sigma**3)*to_pascal):(($10*epsilon/sigma**3)*to_pascal) with yerrorbars linestyle 4 title 'Pressure'

set output "MSD.png"
set title "MSD vs. density in a Lennard-Jones liquid at k_B T=1.2"
set ylabel "MSD"

set xzeroaxis

plot data_file using 1:11:12 with yerrorbars linestyle 4 title 'MSD'

set output "MSD_evol.png"
raw_file_005 = "thermodynamics_0.05.dat"
raw_file_080 = "thermodynamics_0.80.dat"

set xlabel "Time (s)"
set ylabel "MSD (cm^2)"
set title "MSD vs. time in a Lennard-Jones liquid at k_B T=1.2"

f(x) = a*x + b
g(x) = c*x + d

timestep = 10e-4
m2_to_cm2 = 10e4
sigma2 = (3.4e-10)**2

# Fit the data from raw_file_005 to the function f(x)
fit f(x) raw_file_005 using ($0*1000*timestep):($6*sigma2*m2_to_cm2) via a, b

# Fit the data from raw_file_080 to the function g(x)
fit g(x) raw_file_080 using ($0*1000*timestep):($6*sigma2*m2_to_cm2) via c, d

plot raw_file_005 u ($0*1000*timestep):($6*sigma2*m2_to_cm2) t "rho = 0.05", f(x) t'', raw_file_080 u ($0*1000*timestep):($6*sigma2*m2_to_cm2) t "rho = 0.80", g(x) t''

print a/6
print c/6


