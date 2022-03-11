from math import *
from Modules.Donnee import *
import Modules.Cinematique as cm

#--------------------------------------------------------------------------------------------
def Bernstein(m, i, x):
    """ Renvoie le polynome de bernstein (m,i) en x """

    return (factorial(m)/(factorial(i)*factorial(m-i)))*(x**i)*((1-x)**(m-i))

#--------------------------------------------------------------------------------------------
def bezier_curve_n(pts, l):
    n = len(pts)
    res = []
    def B(t):
        P = 0
        for i in range(n):
            P += Bernstein(n, i, t)*pts[i]
        return P
    t = 0
    while t <= l:
        res.append(B(t))
        t += 1/l
    return res

#--------------------------------------------------------------------------------------------
def bezier_point(pts, t):

    if len(pts) == 1:
        return pts[0]
    control_linestring = zip(pts[:-1], pts[1:])
    return bezier_point([(1-t)*p1+t*p2 for p1, p2 in control_linestring],t)

#--------------------------------------------------------------------------------------------
def bezier_curve(pts, nb_pts, l): # algo de Casteljau
    """ Courbe de Bezier pour les n+1 points de controles P, l determine le lissage """

    return [bezier_point(pts,t) for t in (i/((nb_pts-1)*l) for i in range((nb_pts)*l))]

#--------------------------------------------------------------------------------------------
# def h_path(X, Z, Lx, Lz):
#     """ Choisis la bonne hauteur de la courbe de Bezier pour avoir un chemin tjr au dessus du relief """

#     Px, Pz = PtControle(X, Z, Lx, Lz)
#     Z =  bezier_curve(Pz, len(Pz), l)
#     absZ = [abs(k) for k in Z]
#     absPz = [abs(k) for k in Pz]
#     j = absPz.index(max(absPz))
#     for i in range(len(Pz)):
#         if min(absZ) + abs(Pz[i]) > max(absPz) and Pz[i] < Pz[j]:
#             j = i
#     return j

#--------------------------------------------------------------------------------------------
def pt_controle(X0, Z0, Lx, Lz):
    """ On definit les points de controles de la courbe de Bezier par:
    la position initiale, le point le plus haut sur le chemin, et le milieu de la zone d'atterissage.
    Idée: Si on aborde un aplomb on trilatere pour obtenir un point de passage interessant """

    Lz = [-k for k in Lz] # on remet Lz en >= 0
    imax = 0
    tabx, tabz = [X0], [Z0]
    X, Z = [], []

    for i in range(len(Lz)-1): # on determine la zone d'aterissage et le point final
        if Lz[i] == Lz[i+1]:
            xf = (Lx[i]+Lx[i+1])/2
            zf = Lz[i]
            zone = [Lx[i], Lx[i+1]]

    # for i in range(len(Lz)): # si il y a un aplomb on recherche les points concernés et on trilatere
    #     if (zone[0] <= Lx[i] <= zone[1]) and Lz[i] != zf:
    #         X.append(Lx[i])
    #         Z.append(Lz[i])
    # if X:
    #     x, y = trilateration(X,Z)
    #     tabx.append(x)
    #     tabz.append(z)

    for i in range(len(Lz)): # on recupere le plus haut pic sur le chemin en considerant le sens du chemin
        if X0 <= xf and (X0 <= Lx[i] <= xf) and Lz[i] > Lz[imax]:
            imax = i
        elif X0 >= xf and (xf <= Lx[i] <= X0) and Lz[i] > Lz[imax]:
            imax = i
    if Lz[imax] != zf: # si il n'y a pas de relief génant sur le chemin on n'ajoute rien
        tabx.append(Lx[imax])
        tabz.append(Z0) #Lz[imax]
    tabx.append(xf)
    tabz.append(zf)
    #print(tabx, tabz)
    #return Lx, Lz
    return tabx, tabz


#--------------------------------------------------------------------------------------------
def atterissage(power, v_speed, h_speed):
    """ Deuxieme phase à activer une fois au dessus de la zone """

    if abs(v_speed) >= h_speed_max or abs(h_speed) > v_speed_max:
        while power < p_max:
            power += p_var
    else:
         while power > p_min:
            power -= p_var
    print(0, power)

#--------------------------------------------------------------------------------------------
def asservissement(B, m): # ne marche pas du tout
    """ Relie les points de la courbe de Bezier par des segments et renvoie l'angle associé (Ox) """

    dist = (h**2 + w**2)**(1/2) # distance max
    power = 0
    rotate = 0
    X, Z, Xs, Zs = m[0], m[1], m[2], m[3]
    n = len(B)
    for k in range(1,n):
        if B[k-1][0] <= X <= B[k][0]:
            p = (B[k-1][1]-B[k][1])/(B[k-1][0]-B[k][0]) # pente d'un segment de la courbe de bezier
            rotate = pi/2-atan(p) # l'angle est dirigé par rapport a la normale au sol
            print(p)
            for p in range(p_max+1):
                X, Z = cm.Eq_Horaires(p, rotate, Xs, Zs, X, Z, dt)
                Xs, Zs = cm.Vitesse(power, rotate, Xs, Zs, dt)
                m = [round(X, a), round(Z, a), round(Xs, a), round(Zs, a)]
                d = ((X-B[k][0])**2+(Z-B[k][1])**2)**(1/2)
                #d = ((X-(zone[0]+zone[1])/2)**2+(Z-H)**2)**(1/2)
                #print(d)
                if dist >= d:
                    dist = d
                    power = p
    return power, round(rotate, a), m

#--------------------------------------------------------------------------------------------
def trajectoire_approx(B, m):
    X = m[0]
    res = []
    i = 0
    while not zone[0] <= X <= zone[1] and i<80:
        p, r, m = asservissement(B, m)
        X = m[0]
        res.append((p,r))
        i += 1
    #print(res)
    return res






