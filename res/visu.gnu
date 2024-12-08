set terminal pngcairo size 800,600
set output '../res/graph_vitesse.png'
set xlabel "Temps (t)"
set ylabel "Vitesse (v) [m/s]"
set logscale x  # Echelle logarithmique pour l'axe x (temps)
plot "../res/spray.dat" using 1:3 with lines title "v(t) [m/s]"

set output

set terminal png size 800,600
set output '../res/temperature.png'
set xlabel "Temps (t)"
set ylabel "Rayon (r) [µm]"
plot "../res/spray.dat" using 1:6 with lines title "T(t) [°C]"

set output  

set terminal pngcairo size 800,600
set output '../res/graph_rayon.png'
set xlabel "Temps (t)"
set ylabel "Rayon (r) [µm]"
plot "../res/spray.dat" using 1:4 with lines title "r(t) [µm]"

set output  

set terminal pngcairo size 800,600
set output '../res/graph_masse.png'
set xlabel "Temps (t)"
set ylabel "Masse (m) [1e-9 kg]"
plot "../res/spray.dat" using 1:5 with lines title "m(t) [kg]"

set output

set terminal pngcairo size 800,600
set output '../res/humidite.png'
set xlabel "Temps (t)"
set ylabel "Humidité relative [%]"
plot "../res/spray.dat" using 1:8 with lines title "Humidité relative (t)"

set output

set terminal pngcairo size 800,600
set output '../res/vitesse_rel.png'
set xlabel "Temps (t)"
set ylabel "U_s + U_{air} [m/s]"
plot "../res/spray.dat" using 1:9 with lines title "Vitesse relative (t)"

set output

set terminal pngcairo size 800,600
set output '../res/chaleur_sensible.png'
set xlabel "Temps (t)"
set ylabel "T_s + T_{air} [°C]"
plot "../res/spray.dat" using 1:10 with lines title "Chaleur sensible (t)"

set output

