from math import *
import random as rdm

from Modules.Donnee import *
import Modules.Cinematique as cm

#===========================================================================================================#

# Initialisation #--------------------------------------------------------------------------

def rdmpower(power): # O(1)
    """ Renvoie une nouvelle pousée differente de au plus 1 m.s-2 avec l'ancienne """

    if power > p_min and power < p_max:
        return power + rdm.randint(-p_var, p_var)
    elif power == p_max:
        return power + rdm.randint(-p_var, 0)
    else:
       return power + rdm.randint(0, p_var)

#-------------------------------------------------------------------------------------------
def rdmrotate(rotate): # 0(1)
    """ Renvoie un angle different de au plus π/12 avec l'ancien """

    if rotate >= radians(-(r_max-r_var)) and rotate <= radians(r_max-r_var):
        return rotate + radians(rdm.randint(-r_var, r_var))
    elif rotate > radians(r_max-r_var) and rotate <= radians(r_max):
        return rotate + radians(rdm.randint(-r_var, 0))
    else:
       return rotate + radians(rdm.randint(0, r_var))

#-------------------------------------------------------------------------------------------
def creer_individu(p0, r0): # 0(n)
    """ Genere un individu qui respecte les contraintes à partir d'une condition initiale (a t=0) """

    individu = [[p0, r0]]*nb_gene # le premier element de la liste est l'indice d'atterrissage initialisé à 0
    for k in range(nb_gene-1):
        individu[k+1] = [rdmpower(individu[k][0]), rdmrotate(individu[k][1])]
    return individu

#-------------------------------------------------------------------------------------------
def creer_population(p0, r0): # O(n*m)

    pop = [None]*taille_population
    for k in range(taille_population):
        pop[k] = creer_individu(p0, r0)
    return pop

# Mutations #-------------------------------------------------------------------------------

def mutation_individu(ind): # O(n)
    """ Modifie aleatoirement un gene en respectant les contraintes du mars lander """

    if ind == []:
        return []

    n = len(ind)
    for i in range(1, n-2): # indice des élements modifiés: on ne modifie pas la CI
        if rdm.random() < proba_de_muter:
            p, r = ind[i][0], degrees(ind[i][1])
            p1, p2 = ind[i-1][0], ind[i+1][0]
            r1, r2 = degrees(ind[i-1][1]), degrees(ind[i+1][1])
            pl, rl = [p], [r] # on garde intentionellement la poussée non muté car il y a des cas ou la mutation est impossible (ex: [2,3,4])
            for k in range(-p_max, p_max+1):
                if (p1-1 <= p+k <= p1+1) and (p2-1 <= p+k <= p2+1) and p+k != r:
                    pl.append(p+k)
            for k in range(-r_max, r_max+1):
                if (r1-1 <= r+k <= r1+1) and (r2-1 <= r+k <= r2+1):
                    rl.append(r+k)
            ind[i] = [rdm.choice(pl), radians(rdm.choice(rl))]

#-------------------------------------------------------------------------------------------
def mutation_population(pop): # O((n+C(traj))*m)

    for k in range(1, taille_population):
        mutation_individu(pop[k])
        i = cm.Trajectoire(m0, pop[k], dt, fuel, Lx, Lz)[4]
        if abs(degrees(pop[k][i][1])) <= r_var: # on fait muter dès que possible l'angle final du lander afin d'obtenir une solution plus rapidement
            pop[k][i][1] = 0

# Croisement #------------------------------------------------------------------------------

def croiser(ind1, ind2): # O(n) (complexité spatiale plus importante que crossover: O(n^2)?)
    """ Croise deux individus en assurant une continuité, on essaie d'abord de couper en deux à l'indice i,
    puis on regarde si on peut avec i-1 ou i+1, etc (methode discrete)"""

    n = len(ind1)
    i = n//2
    croisements_possibles = []
    fils = []

    def critere_croisement(i, ind1, ind2):
        return (abs(ind1[i][0]-ind2[i][0]) <= 1) and (abs(ind1[i][0]-ind2[i][0]) <= radians(15))

    for k in range(n//2):
        if critere_croisement(i-k, ind1, ind2):
            croisements_possibles.append(i-k)
            fils = ind1[:(i-k)] + ind2[(i-k):]
        if critere_croisement(i+k, ind1, ind2):
            croisements_possibles.append(i-k)
            fils = ind1[:(i+k)] + ind2[(i+k):]

    nb_croisements_possibles = len(croisements_possibles)
    if nb_croisements_possibles > 1:
        fils = ind1
        nb_croisements = rdm.randint(1, nb_croisements_possibles)

        for j in range(1, nb_croisements-1):
            fils = fils[:j-1] + ind2[j-1:j] + fils[j:]
    return fils

#-------------------------------------------------------------------------------------------
def crossover(ind1, ind2): # O(n)
    """ Methode de croisement avec float on prend le barycentre de deux individus (methode continue) """

    b = rdm.random()
    n = len(ind1)
    fils1 = [None]*n
    fils2 = [None]*n
    for i in range(n):
        fils1[i] = [round(b*ind1[i][0]+(1-b)*ind2[i][0]), b*ind1[i][1]+(1-b)*ind2[i][1]]
        fils2[i] = [round((1-b)*ind1[i][0]+b*ind2[i][0]), (1-b)*ind1[i][1]+b*ind2[i][1]]
    return fils1, fils2

#===========================================================================================================#
