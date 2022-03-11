from math import *
from Modules.Donnee import *

#===========================================================================================================#

#-------------------------------------------------------------------------------------------
def Vitesse(power, rotate, Xspeed, Zspeed, dt): # O(1)
    """ Simple pfd + projeté (power vu comme une acceleration initiale """

    Xspeed += -1 * (power * sin(rotate)) * dt
    Zspeed += (power * cos(rotate) - g) * dt
    return Xspeed, Zspeed

#-------------------------------------------------------------------------------------------
def Eq_Horaires(power, rotate, Xspeed, Zspeed, X, Z, dt): # O(1)
    """ On integre pour obtenir les équations horaires (constantes d'intégrations en parametres) """

    X += -1/2 * (power * sin(rotate)) * (dt**2) - Xspeed * dt
    Z += 1/2 * (power * cos(rotate) - g) * (dt**2) + Zspeed * dt
    return X, Z # On arrondit les valeur de x et z

#-------------------------------------------------------------------------------------------
def Trajectoire(m0, ind, dt, fuel, Lx, Lz): # O(n+C(intersection))
    """ Renvoie un tableau constitué des X et Z issu d'un individu (liste de couples (p,r) cf genetique) ainsi que ses parametres d'atterrissage avec m = [X, Z, Xspeed, Zspeed] """
    Lx, Lz = sol_securise(Lx, Lz)
    Lz = [-k for k in Lz]
    m = [k for k in m0] # copie profonde
    X, Z = m[0], m[1]
    Xspeed, Zspeed = m[2], m[3]
    tabX, tabZ = [X], [Z]
    Xp = X
    i = 1
    #b = False
    acc = pts_par_arc # changement de (power, rotate) lorsque acc = pts_par_arc: cela permet de respecter la contrainte changement toutes les secondes et d'avoir plus de point pour la trajectoire

    while (0 <= X <= w) and (0 <= abs(Z) <= h) and (i < len(ind)-1) and fuel >= 0 and not Intersection(Xp, X, Z, Lx, Lz):

        if acc == pts_par_arc: # lisser la trajectoire ==> meilleure approximation de l'arrêt
            power, rotate = ind[i][0], ind[i][1]
            i += 1
            fuel -= power # modelise la consommation en fuel
            acc = 0
        acc += 1
        Xspeed, Zspeed = Vitesse(power, rotate, Xspeed, Zspeed, dt/pts_par_arc)
        X, Z = Eq_Horaires(power, rotate, Xspeed, Zspeed, X, Z, dt/pts_par_arc)
        Xp = tabX[-1]
        # if b:
        #     tabZ.append()
        # else:
        #     tabZ.append(Z)
        tabX.append(X), tabZ.append(Z)
        m = [X, Z, Xspeed, Zspeed]


    return tabX, tabZ, m, fuel, i # i est l'indice d'atterrisage

# fonctions de detections des collisions #----------------------------------------------------

def Intersection(Xp, X, Z, Lx, Lz):
    """ Materialise le sol et donne des conditions d'arrêt quand on a un point en dessous du sol """

    ptc = []
    Lz = [-k for k in Lz]
    def segment(X, i):
        a = ((Lz[i]-Lz[i+1])/(Lx[i]-Lx[i+1])) # on a multiplié par -1 pour l'axe z
        b = (Lz[i]*Lx[i+1]-Lz[i+1]*Lx[i])/(Lx[i+1]-Lx[i])
        return (a*X+b, a, b)

    for i in range(len(Lx)-1):
        if (Lx[i] <= X < Lx[i+1]):
            ptc.append(i)

    if len(ptc) == 1:
        return abs(Z) <= segment(X, ptc[0])[0]
    elif len(ptc) == 2:
        i,j = ptc
        for k in ptc:
            if min(Lz[k], Lz[k+1]) <= abs(Z) <= max(Lz[k], Lz[k+1]) and abs(segment(X, k)[1]) > 1:
                z, a, b = segment(X , k)
                print("coucou", k)
                return (Xp <= (z-b)/a <= X) or (X <= (z-b)/a <= Xp) and abs(Z) <= z
        if max(Lz[i], Lz[i+1]) <= min(Lz[j], Lz[j+1]):
            # if abs(Z) <= max(Lz[i], Lz[i+1]) and abs(segment(X, i)[1]) > 1:
            #     z, a, b = segment(X , i)
            #     return (z-b)/a
            # elif abs(Z) >= max(Lz[i], Lz[i+1]) and abs(segment(X, i)[1]) > 1:
            #     return
            #else:
            #print(not (segment(X, i)[0] <= abs(Z) <= segment(X, j)[0]))
            #print(not (segment(X, i)[0] <= abs(Z) <= segment(X, j)[0])) or (not (segment(X, j)[0] <= abs(Z) <= segment(X, i)[0]))
            return not (segment(X, i)[0] <= abs(Z) <= segment(X, j)[0])
        elif max(Lz[j], Lz[j+1]) <= min(Lz[i], Lz[i+1]):
            #print(not (segment(X, j)[0] <= abs(Z) <= segment(X, i)[0]))
            #print(not (segment(X, i)[0] <= abs(Z) <= segment(X, j)[0])) or (not (segment(X, j)[0] <= abs(Z) <= segment(X, i)[0]))
            return not (segment(X, j)[0] <= abs(Z) <= segment(X, i)[0])

        #elif min(Lz[j], Lz[j+1]) <= abs(Z) <= max(Lz[j], Lz[j+1]) and abs(segment(X, j)[1]) > 1:

        return (not (segment(X, i)[0] <= abs(Z) <= segment(X, j)[0])) or (not (segment(X, j)[0] <= abs(Z) <= segment(X, i)[0]))

        # else:
        #     if X != Lx[i+1]:
        #         print("sol non valide")

    else:
        print("zone non accessible")
        #return True
        #print ("Lander non détecté")

#----------------------------------------------------------------------------------------------
def collision(Xspeed, Zspeed, rotate): # O(1)
    """ Conditions d'atterissages du lander """

    return (round(degrees(rotate)) == 0) and abs(Xspeed) <= v_speed_max and abs(Zspeed) <= h_speed_max




