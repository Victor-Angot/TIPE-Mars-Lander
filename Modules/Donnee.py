from math import *

# Parametres graphiques #-----------------------------------------------------------------------

scale = 0.1
w, h = 7000 * scale, 3000 * scale  # hauteur et largeur de la fenetre

# Constantes et conditions initiales #----------------------------------------------------------

g = 3.711 # pesanteur martienne

X0 = 2500 * scale
Z0 = 2700 * scale

Xspeed0 = 0
Zspeed0 = 0

fuel = 550

rotate0 = 0
power0 = 0

m0 = [X0 , Z0, Xspeed0, Zspeed0]

# Contraintes #---------------------------------------------------------------------------------

r_max, r_min = 90, -90 # angle maximal (en degré)
r_var = 15 # rotation maximale d'un tour à l'autre (en degré)
p_max, p_min = 4, 0 # poussée maximale (en m/s^2)
p_var = 1
dt = 1 # pas de temps (ici changement de (p,r) toutes les secondes)
h_speed_max, v_speed_max = 40, 20 # en m/s

# Parametres AG #--------------------------------------------------------------------------------

proba_de_muter = 0.05 # probabilité d'un individu (chromosome) de muter
nb_gene = 40 # taille d'un individu (notée n)
taille_population = 80 # nombre d'individu d'une population (notée m)
nb_de_generation_max = 1000 # notée N

pourcentage_grade_retenu = 0.2 # param de la selection elitiste
nb_grade_retenu = int(taille_population * pourcentage_grade_retenu)
proba_non_grade_retenu = 0.05

proba_de_selection = 0.8 # proba de l'indivu à la fitness forte de gagner dans la selection par tournois
nb_parents = taille_population * 0.4

# Autres #--------------------------------------------------------------------------------------

graph = []
note_moyenne = []
a = 3  # precision de l'arrondi des valeurs numériques
l = 10 # lissage des courbes de Bezier
pts_par_arc = 5 # lissage des trajectoires (10 lent mais très bon resultat)
delta = 0 # rehaussage du sol pour des raisons de sécurité

def sol_securise(Lx,Lz):
    return Lx, [delta - k for k in Lz]

# Base de données: generation des sols #---------------------------------------------------------

def landing_zone(Lx, Lz):
    """ Renvoie une liste contenant les deux extremités de la zone d'atterissage """
    for i in range(len(Lz)-1):
        if Lz[i] == Lz[i+1]:
            return [Lx[i], Lx[i+1]]
    print("Pas d'atterissage possible")


#-------------------------------------------------------------------------------------------
def h_zone(Lx, Lz):
    """ Renvoie la hauteur z de la zone d'atterissage """

    for i in range(len(Lz)-1):
        if Lz[i] == Lz[i+1]:
            return abs(Lz[i])


Coord_X = [[0, 1000, 1500, 3000, 4000, 5500, 6999],
[0, 1000, 1500, 3000, 3500, 3700, 5000, 5800, 6000, 6999],
[0, 1000, 1500, 3000, 4000, 5500, 6999],
[0, 300, 350, 500, 800, 1000, 1200, 1500, 2000, 2200, 2500, 2900, 3000, 3200, 3500, 3800, 4000, 5000, 5500, 6999],
[0, 300, 350, 500, 1500, 2000, 2500, 2900, 3000, 3200, 3500, 3800, 4000, 4200, 4800, 5000, 5500, 6000, 6500, 6999],
[0, 300, 1000, 1500, 1800, 2000, 2200, 2400, 3100, 3150, 2500, 2200, 2100, 2200, 3200, 3500, 4000, 4500, 5000, 5500, 6000, 6999],
[0, 300, 1000, 2000, 2500, 3700, 4700, 4750, 4700, 4000, 3700, 3750, 4000, 4900, 5100, 5500, 6200, 6999]]
# Ensemble des Lx

Coord_Z = [[100, 500, 1500, 1000, 150, 150, 800],
[100, 500, 100, 100, 500, 200, 1500, 300, 1000, 2000],
[100, 500, 1500, 1000, 150, 150, 800],
[1000, 1500, 1400, 2000, 1800, 2500, 2100, 2400, 1000, 500, 100, 800, 500, 1000, 2000, 800, 200, 200, 1500, 2800],
[1000, 1500, 1400, 2100, 2100, 200, 500, 300, 200, 1000, 500, 800, 200, 800, 600, 1200, 900, 500, 300, 500],
[450, 750, 450, 650, 850, 1950, 1850, 2000, 1800, 1550, 1600, 1550, 750, 150, 150, 450, 950, 1450, 1550, 1500, 950, 1750],
[1800, 1200, 1550, 1200, 1650, 220, 220, 1000, 1650, 1700, 1600, 1900, 2100, 2050, 1000, 500, 800, 600]]
# Ensemble des Lz

Ci = [[250, 270, 0, 0], [650, 280, 100, 0], [650, 280, 90, 0], [50, 270, -100, 0], [650, 270, 50, 0], [600, 260, 20, 0], [650, 200, 0, 0]] # ensemble des [X0 , Z0, Xspeed0, Zspeed0]

carburant = [550, 600, 750, 800, 1000, 1000, 1200]

gene0 = [(0, 0), (0, radians(-90)), (0, radians(-90)), (0, radians(90)), (0, radians(-90)), (0, radians(-45)), (0, 0)] # (p,r) initiale


# Parametres initiaux #------------------------------------------------------------------------------------

Lx = [k * scale for k in Coord_X[0]]
Lz = [-k * scale for k in Coord_Z[0]]

zone = landing_zone(Lx, Lz)
H = h_zone(Lx, Lz)
largeur_zone = zone[1] - zone [0]
distance_milieu_zone = (zone[1] + zone[0])/2 - zone[0]
lnote_moyenne = []
lnote_max = []

