# Premier graphique : vitesse et température
set terminal pngcairo size 800,600
set output '../res/graph_vitesse_temperature.png'

# Configuration des axes
set xlabel "Temps (t)"
set ylabel "Vitesse (v) [m/s]"
set logscale x  # Echelle logarithmique pour l'axe x (temps)

# Définir l'axe secondaire pour la température
set y2label "Température (T) [°C]"
set y2range [14:21]  # Plage de température
set y2tics 1  # Espacement des tics pour la température

# Tracer les courbes de vitesse et température
plot "../res/drop_0.dat" using 1:3 with lines title "v(t) [m/s]", "../res/drop_0.dat" using 1:6 axes x1y2 with lines title "T(t) [°C]"

# Réinitialiser la sortie PNG pour le second graphique
set output  # Annuler la sortie PNG en cours

# Deuxième graphique : rayon et masse
set terminal pngcairo size 800,600
set output '../res/graph_rayon_masse.png'

# Configuration des axes
set xlabel "Temps (t)"
set ylabel "Rayon (r) [µm]"

set yrange [40:110]

# Réinitialiser les axes secondaires pour le second graphique
set y2label "Masse (m) [1e-9 kg]"
set y2range [0:5]  
set y2tics 0.5  

# Tracer les courbes de rayon et masse
plot "../res/drop_0.dat" using 1:4 with lines title "r(t) [µm]", "../res/drop_0.dat" using 1:5 axes x1y2 with lines title "m(t) [kg]"

# Réinitialiser la sortie PNG en cours (non nécessaire ici, mais pour la propreté du code)
set output
