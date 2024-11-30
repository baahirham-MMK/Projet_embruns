# Compilateur utilisé
CC = nvc++

# Options en mode optimisé - La variable DEBUG est définie comme fausse
OPTIM_FLAG = -O3 -DNDEBUG -std=c++11 -acc -gpu=managed -Minfo=accel
# Options en mode debug - La variable DEBUG est définie comme vraie
DEBUG_FLAG = -g -DDEBUG -std=c++11 -acc -gpu=managed -Minfo=accel

# On choisit comment on compile
CXX_FLAGS = $(DEBUG_FLAG)

# Le nom de l'exécutable
PROG = run

# Les fichiers source à compiler
SRC = main.cpp DataFile.cpp Function.cpp Spray_openacc.cpp

# La commande complète : compile seulement si un fichier a été modifié
$(PROG): $(SRC)
	$(CC) $(SRC) $(CXX_FLAGS) -o $(PROG)

# Évite de devoir connaître le nom de l'exécutable
all: $(PROG)

# Supprime l'exécutable, les fichiers binaires (.o) et les fichiers temporaires de sauvegarde (~)
clean:
	rm -f *.o ../res/*.png ../res/*.dat *~ $(PROG)

