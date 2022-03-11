# Librairies #-------------------------------------------------------------------------------

from math import *
from tkinter import *
from tkinter import messagebox
import numpy as np
from time import *
import matplotlib.pyplot as plt
import random as rdm

from Modules.Donnee import *
import Modules.Bezier as bz
import Modules.Cinematique as cm
import Modules.Genetique as gt

#===========================================================================================================#


# Fonctions utilitaires #--------------------------------------------------------------------

def zip(L1, L2):
    """ transforme de liste X, Z en une liste contenant les couples (x, -z) """

    tab = []
    for i in range(len(L1)):
        tab.append((L1[i], -L2[i])) # signe moins pour representation graphique (l'axe z est orienté vers le bas)
    return tab

#-------------------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------------------
def distance(ind1, ind2):

    X1, Z1 = cm.Trajectoire(m0, ind1, dt, fuel, Lx, Lz)[0], cm.Trajectoire(m0, ind1, dt, fuel, Lx, Lz)[1]
    X2, Z2 = cm.Trajectoire(m0, ind2, dt, fuel, Lx, Lz)[0], cm.Trajectoire(m0, ind2, dt, fuel, Lx, Lz)[1]
    distance_collision = ((X1[-1]-X2[-1])**2+(Z1[-1]-Z2[-1])**2)**(1/2)
    return distance_collision

def repart_distance(n):
    tab = [0] * 700
    for k in range(n):
        ind1 = gt.creer_individu(0,0)
        ind2 = gt.creer_individu(0,0)
        tab[round(distance(ind1, ind2))] += 1
    X = [k for k in range(700)]
    tab = [k/n for k in tab]
    plt.bar(X, tab, align="center")
    plt.show()

#-------------------------------------------------------------------------------------------
def distance_au_sol(X, Z, Lx, Lz):
    """ Calcule la longueur des segments du sol cumulés de X à la zone d'atterissage """

    Lz = [-k for k in Lz]
    p, q = 0, 0 # indice du segment considéré

    def dist_eucli(a, b):
        return ((a[0]-b[0])**2+(a[1]-b[1])**2)**(1/2)

    mid_zone = (zone[0] + zone[1])/2
    x = (mid_zone, H)
    for i in range(len(Lx)-1):
        if (Lx[i] <= X <= Lx[i+1]):
            p = i
        if (Lx[i] <= mid_zone <= Lx[i+1]):
            q = i

    if p <= q:
        d = dist_eucli((X,Z), (Lx[p+1], Lz[p+1])) # distance de X au premier point de discontinuité
        for i in range(p+1, q+1):
            a = (Lx[i], Lz[i])
            b = (Lx[i+1], Lz[i+1])
            d += dist_eucli(a, b)
    else:
        d = dist_eucli((X,Z), (Lx[p-1], Lz[p-1]))
        for i in range(p-1, q, -1):
            a = (Lx[i], Lz[i])
            b = (Lx[i+1], Lz[i+1])
            d += dist_eucli(b, a)
    return d + dist_eucli((Lx[q], Lz[q]), x)

#===========================================================================================================#

# Parametres #-----------------------------------------------------------------------------

proba_de_muter = 0.05 # probabilité d'un individu (chromosome) de muter
nb_gene = 40 # taille d'un individu (que l'on notera n pour la complexité)
taille_population = 80 # nombre d'individu d'une population
nb_de_generation_max = 150

pourcentage_grade_retenu = 0.2 # param de la selection elitiste
nb_grade_retenu = int(taille_population * pourcentage_grade_retenu)
proba_non_grade_retenu = 0.05

proba_de_selection = 0.6 # proba de l'indivu à la fitness forte de gagner dans la selection par tournois
nb_parents = taille_population * 0.4

# Evaluation #------------------------------------------------------------------------------

def scaling(note, k):
    """ scaling exponentielle """

    return note**k

#-------------------------------------------------------------------------------------------
def fs(n, p):
    """ n génération de la population, p paramétre, cette fonction donne k"""

    return (tan((n/(nb_de_generation_max+1))*(pi/2)))**p

#-------------------------------------------------------------------------------------------
def S(d, sig, a):
    """ sig voisinage de Xc, a un parametre """
    if d < sig:
        return 1-(d/sig)**a
    else:
        return 0.01


def sharing(note, individu, population):
    """ sharing non optimisé (on calcule la distance à chaque individu) """

    # ameliorer en calculant le Xc une seule fois

    res = 0
    for k in population:
       res += S(distance(k, individu), 50, 2)
    note = note/res
    return note

#-------------------------------------------------------------------------------------------
def evalue(ind): # C(traj)
    """ Attribue un score proportionelle à la distance à la zone d'atterissage et
    l'angle/vitesse d'arrivé(e) que l'on cherche à minimiser. Renvoie un boléen selon si l'individu
    est éligible à être solution """

    X, Z, m, fuel_f, i = cm.Trajectoire(m0, ind, dt, fuel, Lx, Lz)
    Xc, Zc = X[-1], Z[-1] # point après collision/atterrissage
    Xspeed, Zspeed = m[2], m[3]
    mid_zone = (zone[0]+zone[1])/2
    L = zone[1]-zone[0]
    rf = ind[i][1]

    note_distance = 100*((abs(Xc-mid_zone))/(abs(X0-mid_zone)))
    # rapport distance au mileu de la zone d'atterrissage par rapport sur distance initiale

    #note_distance = distance_au_sol(Xc, Zc, Lx, Lz)
    note = 500
    if (zone[0]+L*0.05 <= Xc <= zone[1]-L*0.05): #pénalités
        #note_distance = max(note_distance/10, 50) # on reduit l'importance de la note_distance
        if abs(Xspeed) >= v_speed_max:
            note -= 100*(abs(Xspeed)/300) # on pose une vitesse maximale sinon le modéle perd son sens
        if abs(Zspeed) >= h_speed_max:
            note -= 100*(abs(Zspeed)/300)
        if abs(degrees(rf)) >= r_var:
            note -= 100*(abs(round(degrees(rf)))/90)
    note -= note_distance
    note /= 5

    est_sol = (zone[0]+L*0.05 <= Xc <= zone[1]-L*0.05) and cm.collision(Xspeed, Zspeed, rf) and Z[-1] <= H <= Z[-2] # on ressere la zone d'atterrissage de maniere arbitraire pour éviter les cas limites, génants en pratique
    return note, est_sol #scaling(note, fs(int(cpt.get()), 1))

#-------------------------------------------------------------------------------------------
def evalue_population(population):

    ind_note = [None]*taille_population
    k = 0
    for ind in population:
        note, est_sol = evalue(ind)
        ind_note[k] = (note, ind, est_sol)
        k += 1
    return sorted(ind_note, reverse = True)

# Selection et évolution #-------------------------------------------------------------------

def evolution_opti(pop):
    npop = evalue_population(pop)
    solution = []
    pop_evaluee = [None]*taille_population
    k = 0
    for note, ind, est_sol in npop:
        if est_sol:
            solution.append(ind)
        pop_evaluee[k] = ind
        k += 1
    if solution != []:
        return [pop_evaluee, solution]
    parents = pop_evaluee[:nb_grade_retenu]
    for ind in pop_evaluee[nb_grade_retenu:]:
        if rdm.random() < proba_non_grade_retenu:
            parents.append(ind)
    while len(parents) < nb_parents:
        ind1, ind2 = rdm.choice(npop), rdm.choice(npop)
        if ind1[0] > ind2[0]:
            if rdm.random() < proba_de_selection:
                parents.append(ind1[1])
        else:
            if rdm.random() < proba_de_selection:
                parents.append(ind2[1])
    fils = []
    n = nb_gene
    while len(fils) < taille_population-len(parents):

        ind1, ind2 = rdm.choice(parents), rdm.choice(parents)
        b = rdm.random()
        fils1 = [None] * n
        fils2 = [None] * n
        for i in range(n):
            fils1[i] = [round(b*ind1[i][0]+(1-b)*ind2[i][0]), b*ind1[i][1]+(1-b)*ind2[i][1]]
            fils2[i] = [round((1-b)*ind1[i][0]+b*ind2[i][0]), (1-b)*ind1[i][1]+b*ind2[i][1]]
        fils.append(fils1)
        fils.append(fils2)
    fils = parents + fils

    return [fils, solution]

#--------------------------------------------------------------------------------------------

def evolution(pop):
    """ Passage de la population n à n+1 """

    npop = evalue_population(pop)
    somme_note = 0
    solution = []
    pop_evaluee = [None]*taille_population
    k = 0

    # On separe les notes et les individus, de plus on donne un critere de convergence

    for note, ind, est_sol in npop:
        #note = sharing(note, ind, pop[0])
        somme_note += note
        pop_evaluee[k] = ind
        if est_sol: # Le lander est sur la zone d'atterrisage et respecte les contraintes
            solution.append(ind)
        k += 1

    note_moyenne = somme_note/taille_population
    note_max = npop[0][0]

    if solution:
        return pop_evaluee, [note_moyenne, note_max], solution

    """ Methode élitiste: """
    # On conserve les individus les mieux notés et un nombre aleatoire d'individus plus faibles pour eviter de      converger vers un maximum local

    parents = pop_evaluee[:nb_grade_retenu]
    for ind in pop_evaluee[nb_grade_retenu:]:
        if rdm.random() < proba_non_grade_retenu:
            parents.append(ind)

    """ Selection par tournois """

    # #parents = []
    # while len(parents) < nb_parents:
    #     ind1, ind2 = rdm.choice(npop), rdm.choice(npop)
    #     if ind1[0] > ind2[0]: # on cherche à maximiser la note
    #         if rdm.random() < proba_de_selection:
    #             parents.append(ind1[1])
    #     else:
    #         if rdm.random() < proba_de_selection:
    #             parents.append(ind2[1])

    """ Roulette wheel """ # O(nlog(n))

    roulette = [None]*taille_population

    k = 0
    for note, ind, est_sol in npop:
        roulette[k] = [(note/somme_note), k]
        k += 1
    roulette = sorted(roulette, reverse = True) # on trie par ordre décroissant
    for i in range(len(roulette)):
        for j in range(i+1, len(roulette)):
            roulette[i][0] += roulette[j][0]
    while len(parents) < nb_parents:
        alea = rdm.random()
        for k in range(1, len(roulette)):
            if roulette[k][0] <= alea <= roulette[k-1][0]:
                parents.append(pop_evaluee[roulette[k][1]])

    """ Si la note moyenne ne varie pas (extrema local) on conserve les meilleurs individus et on génére une nouvelle population pour completer les parents """

    n = len(lnote_moyenne)
    it = 0
    if n >= 5 and sum(lnote_moyenne[n-6:])/(lnote_moyenne[-1]*5) >= 0.95: # si les 5 dernieres note moyenne ne varie pas de plus de 5% (arbitraire)
        parents = parents[:nb_grade_retenu//2]
        new_pop = gt.creer_population(power0, rotate0)
        while len(parents) < nb_parents:
            parents.append(new_pop[it])
            it += 1

    fils = []
    while len(fils) < taille_population-len(parents):
        """ croisement continu """

        ind1, ind2 = rdm.choice(parents), rdm.choice(parents)
        fils1, fils2 = gt.crossover(ind1, ind2)
        fils.append(fils1)
        fils.append(fils2)

        """ croisement discret """
        # ind3 = gt.croiser(ind1, ind2)
        # if ind3 != []:
        #     fils.append(ind3)
        # else:
        #     #print("échec croisement")
        #     fils.append(ind1), fils.append(ind2)

    fils = parents + fils

    # On applique les mutations
    gt.mutation_population(fils)

    return fils, [note_moyenne, note_max], solution


# Fonctions graphiques #---------------------------------------------------------------------

def trace_traj():
    """ Represente la trajectoire (Z en fonction de X) """

    global object
    can.delete(object)
    individu_rdm = gt.creer_individu(power0, rotate0)
    X, Z, m, fuel_f, i = cm.Trajectoire(m0, individu_rdm, dt, fuel, Lx, Lz)
    tab = zip(X, Z)
    #print(evalue(individu_rdm)[0])
    object = can.create_line(tab, smooth = "true", width = 1, fill = "#{:x}{:x}{:x}".format(rdm.randint(150, 255), rdm.randint(150, 255), rdm.randint(150, 255)))
    object

#--------------------------------------------------------------------------------------------
def trace_bezier():
    """ Trace la courbe de Bezier associé au point du sol plus """

    Px, Pz = bz.pt_controle(m0[0], m0[1], Lx, Lz)
    X, Z = bz.bezier_curve(Px, len(Px), l), bz.bezier_curve(Pz, len(Pz), l)
    #X, Z = bz.bezier_curve_n(Px, l), bz.bezier_curve_n(Pz, l)
    tab = zip(X, Z)
    can.create_line(tab, smooth = "true", width = 2, fill = "yellow")

    """ Trace la trajectoire approximé autour de la courbe de Bezier """

    # B = zip(X, [-k for k in Z])
    # ind = bz.trajectoire_approx(B, m0)
    # X, Z, m, fuel_f = cm.Trajectoire(m0, ind, dt, fuel, Lx, Lz)
    # tab2 = zip(X, Z)
    # can.create_line(tab2, smooth = "true", width = 1, fill = "white")

#--------------------------------------------------------------------------------------------
freq = [0]*nb_de_generation_max

def trace_evolution(item = None):
    """ Represente graphiquement les générations de populations """
    global freq
    global pop, id_anim, graph
    global lnote_moyenne, lnote_max
    if graph:
        for individu in graph:
            can.delete(individu)
    graph = []
    max_fuel = 0

    if int(cpt.get()) == 0:
        pop = [gt.creer_population(power0, rotate0), [], []]
    else:
        pop = evolution(pop[0])
        lnote_moyenne.append(pop[1][0])
        lnote_max.append(pop[1][1])

    if pop[2] or int(cpt.get()) >= nb_de_generation_max: # si il existe une solution ou si nombre de génération trop grand
        for sol in pop[2]:
            X, Z, m, fuel_f, i = cm.Trajectoire(m0, sol, dt, fuel, Lx, Lz)
            max_fuel = max(max_fuel, fuel_f)
            tab = zip(X, Z)
            g_sol = can.create_line(tab, smooth = "true", width = 3, fill = "green")
        if pop[2]:
            freq[int(cpt.get())] += 1
            print("Atterrisage réussi: ",fuel, "-->", max_fuel)
        else:
            freq[-1] += 1
            print("Echec de l'atterrisage")

    else:
        for ind in pop[0]:
            X, Z, m, fuel_f, i = cm.Trajectoire(m0, ind, dt, fuel, Lx, Lz)
            tab = zip(X, Z)
            graph.append(can.create_line(tab, smooth = "true", width = 1, fill = "white"))

    for individu in graph:
        individu

    cpt.set(int(cpt.get())+1)
    if not pop[2] and int(cpt.get()) < nb_de_generation_max:
        id_anim = app.after(1, trace_evolution, item)

#--------------------------------------------------------------------------------------------
def stop():
    can.after_cancel(id_anim)

#-------------------------------------------------------------------------------------------
ncv = 0

def frequence(item = None):
    global anim
    global ncv, freq
    restart()
    trace_evolution()
    ncv += 1
    print(ncv)
    if ncv < 100:
        anim = app.after(30000, frequence, item)
    else:
        L=[]
        X = [k for k in range(nb_de_generation_max)]
        freq = [k/ncv for k in freq]
        plt.figure("Nombre de géneration avant convergence pour "+str(ncv)+" essais")
        plt.bar(X, freq, align="center")
        for k in range(nb_de_generation_max):
            if k%10 == 0:
                L.append(k)
            else:
                L.append("")

        plt.xticks(X, L)
        plt.show()

#--------------------------------------------------------------------------------------------
def Niveau1():
    """ Genere les segments qui represente le sol du niveau 1"""

    global Coord_X, Coord_Z, Lx, Lz
    global zone, H, m0, fuel, power0, rotate0
    fuel = carburant[0]
    m0 = Ci[0]
    power0, rotate0 = gene0[0]
    Lx, Lz = [k * scale for k in Coord_X[0]], [-k * scale for k in Coord_Z[0]]
    #Lx, Lz = sol_securise([k * scale for k in Coord_X[0]], [-k * scale for k in Coord_Z[0]])
    #Lz = [-k for k in Lz]
    zone = landing_zone(Lx, Lz)
    H = h_zone(Lx, Lz)
    restart()

#--------------------------------------------------------------------------------------------
def Niveau2():
    global Coord_X, Coord_Z, Lx, Lz
    global zone, H, m0, fuel, power0, rotate0
    fuel = carburant[3]
    # m0 = Ci[3]
    # power0, rotate0 = gene0[3]
    Lx, Lz = [k * scale for k in Coord_X[3]], [-k * scale for k in Coord_Z[3]]
    zone = landing_zone(Lx, Lz)
    H = h_zone(Lx, Lz)
    restart()

#--------------------------------------------------------------------------------------------
def Niveau3():
    global Coord_X, Coord_Z, Lx, Lz
    global zone, H, m0, fuel, power0, rotate0
    fuel = carburant[4]
    # m0 = Ci[4]
    # power0, rotate0 = gene0[4]
    Lx, Lz = [k * scale for k in Coord_X[4]], [-k * scale for k in Coord_Z[4]]
    zone = landing_zone(Lx, Lz)
    H = h_zone(Lx, Lz)
    restart()

#--------------------------------------------------------------------------------------------
def Niveau4():
    global Coord_X, Coord_Z, Lx, Lz
    global zone, H, m0, fuel, power0, rotate0
    fuel = carburant[5]
    m0 = Ci[5]
    power0, rotate0 = gene0[5]
    Lx, Lz = [k * scale for k in Coord_X[5]], [-k * scale for k in Coord_Z[5]]
    zone = landing_zone(Lx, Lz)
    H = h_zone(Lx, Lz)
    restart()

#--------------------------------------------------------------------------------------------
def Niveau5():
    global Coord_X, Coord_Z, Lx, Lz
    global zone, H, m0, fuel, power0, rotate0
    fuel = carburant[6]
    m0 = Ci[6]
    power0, rotate0 = gene0[6]
    Lx, Lz = [k * scale for k in Coord_X[6]], [-k * scale for k in Coord_Z[6]]
    zone = landing_zone(Lx, Lz)
    H = h_zone(Lx, Lz)
    restart()

#--------------------------------------------------------------------------------------------
def info():
    print("************************************")
    print("\nLanding zone: ", zone, H)
    print("Note moyenne: ", round(pop[1][0]))
    note_max, ind, est_sol = evalue_population(pop[0])[0]
    print("\nNote max: ", round(note_max))
    tabX, tabZ, m, fuel_f, i = cm.Trajectoire(m0, ind, dt, fuel, Lx, Lz)
    if pop[2]:
        ind = pop[2][0]
        ind = ind[:i+1]
        print("\nNote solution: ", evalue(ind)[0])
        print(ind)
    print("X_speed: ", round(m[2]), "m/s")
    print("Z_speed: ", round(m[3]), "m/s")
    print("Angle_atterissage: ", round(degrees(ind[i][1])),"°")
    print("X_atterisage: ", round(m[0]))
    print("Z_atterisage: ", (round(tabZ[-2]), round(tabZ[-1])), "\n")
    print("************************************")

#--------------------------------------------------------------------------------------------
def erreur():
    """ Graphe de l'evolution du fitnesse score """

    def make_plot(lnote_moyenne, lnote_max):
        gen = range(1, len(lnote_moyenne)+1)
        plt.figure("Evolution de la note au cours des generations")
        ax1 = plt.subplot(2,1,1)
        ax1.plot(gen, lnote_moyenne, label = "note moyenne")
        ax2 = plt.subplot(2,1,2)
        ax2.plot(gen, lnote_max,color = "red", label = "note max")
        ax1.set_ylabel("Note")
        ax2.set_xlabel("Générations")
        ax2.set_ylabel("Note")
        ax1.grid()
        ax1.legend()
        ax2.grid()
        ax2.legend()
        plt.show()
        plt.close()
    make_plot(lnote_moyenne, lnote_max)

#--------------------------------------------------------------------------------------------
def quit():
    print("Fin de la simulation")
    app.destroy()

#--------------------------------------------------------------------------------------------
def restart():
    """ Redémarre la simulation """

    global lnote_moyenne, lnote_max, ncv
    #print("La simulation a été redémarrée")
    cpt.set('0') # on reinitialise le compteur de l'AG
    lnote_moyenne, lnote_max = [], []
    can.delete('all')
    for i in range(len(Lx)-1):
        can.create_line(Lx[i], Lz[i], Lx[i+1], Lz[i+1], width = 1, fill = "red")
    Lxup, Lzup = sol_securise(Lx, Lz)
    Lzup = [-k for k in Lzup]
    for i in range(len(Lxup)-1):
        can.create_line(Lxup[i], Lzup[i], Lxup[i+1], Lzup[i+1], width = 1, fill = "red", dash = (4, 4))

#--------------------------------------------------------------------------------------------
def a_propos():
    mon_message = messagebox.showinfo("Ce TIPE a été réalisé par :", "Victor Angot \n MP* \n Lycée Massena")

#================================================= MAIN =================================================#

app = Tk()
app.title("Project Mars Lander: T5")

can = Canvas(app, scrollregion = (0, -h, 0, -h), bg = "black", height = h, width = w) # on a deplacé l'origine avec scrollregion
can.pack(side = "top")

mon_menu = Menu(app)
app.config(menu = mon_menu)

# Affichage par defaut du niveau 1
Lx = [k * scale for k in Coord_X[0]]
Lz = [-k * scale for k in Coord_Z[0]]
for i in range(len(Lx)-1):
    can.create_line(Lx[i], Lz[i], Lx[i+1], Lz[i+1], width = 1, fill = "red")
Lxup, Lzup = sol_securise(Lx, Lz)
Lzup = [-k for k in Lzup]
for i in range(len(Lxup)-1):
    can.create_line(Lxup[i], Lzup[i], Lxup[i+1], Lzup[i+1], width = 1, fill = "red", dash = (4, 4))


# Menu fichier
fichier = Menu(mon_menu, tearoff = 0)
fichier.add_command(label = "Redémarrer", command = restart)
fichier.add_separator()
fichier.add_command(label = "Quitter", command = quit)
mon_menu.add_cascade(label = "Fichier", menu = fichier)

# Menu cartes
cartes = Menu(mon_menu, tearoff = 0)
cartes.add_command(label = "Niveau 1", command = Niveau1)
cartes.add_separator()
cartes.add_command(label = "Niveau 2", command = Niveau2)
cartes.add_separator()
cartes.add_command(label = "Niveau 3", command = Niveau3)
cartes.add_separator()
cartes.add_command(label = "Niveau 4", command = Niveau4)
cartes.add_separator()
cartes.add_command(label = "Niveau 5", command = Niveau5)

mon_menu.add_cascade(label = "Affichage", menu = cartes)

# A propos
apropos = Menu(mon_menu, tearoff = 0)
mon_menu.add_cascade(label = "A propos", menu = apropos)
apropos.add_command(label = "TIPE", command = a_propos)

# Boutons
# bou1 = Button(app, text = "Trajectoire aléatoire", command = trace_traj)
# bou1.pack()

# bou2 = Button(app, text = "Courbe de Bezier", command = trace_bezier)
# bou2.pack()

bou6 = Button(app, text = "Graphe erreur", command = erreur)
bou6.pack(side = "top")

bou7 = Button(app, text = "Convergence", command = frequence)
bou7.pack(side = "top")

bou4 = Button(app, text = "Information", command = info)
bou4.pack()

# Compteur de generation
cpt = StringVar()
cpt.set('0')

lbl = Label(app, width = '10', textvariable = cpt, font = 'Avenir 15 bold')
lbl.pack(side = "top", padx = 10, pady = 10)

bou3 = Button(app, text = "Simulation", command = trace_evolution)
bou3.pack(side = "top", padx = 0, pady = 0)

bou5 = Button(app, text = "Pause", command = stop)
bou5.pack(side = "top")

app.mainloop()
