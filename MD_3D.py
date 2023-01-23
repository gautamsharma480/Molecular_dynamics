#from vpython import *
import vpython
from vpython import *
import numpy as np
import random
import matplotlib
# MD code in 2D
Natm=125
Nmax=125
Tini = 12.0
dens = 1.0
t1 = 0

x = np.zeros(Nmax)
#print(x.size)
y = np.zeros(Nmax)
z = np.zeros(Nmax)
vx = np.zeros(Nmax)
vy = np.zeros(Nmax)
vz = np.zeros(Nmax)
#print(vy)
fx = np.zeros((Nmax,2)) # new and old forces are needed, so two dimensions
fy = np.zeros((Nmax,2)) # same reasoning.
fz = np.zeros((Nmax,2))
#print(fx.size,fy.size, fx.shape)
L = int(1*np.cbrt(Natm))
#print('length of the box = ', L)

def twelveran():
    s = 0
    for i in range(1,13):
        s += random.random()
        return s/12 - 0.5
#print(twelveran())
def initial_pos_vel():
    i = -1 # to start indexing from 0.
    for ix in range(0,L):
        for iy in range(0,L):
            for iz in range(0,L):
                i = i + 1 # this changes the index of array x[i], y[i] from 0 to Nmax-1.
                #print(ix,iy,iz,i)
                x[i] = ix
                y[i] = iy
                z[i] = iz
                vx[i] = twelveran()
                vy[i] = twelveran()
                vz[i] = twelveran()
                vx[i] = vx[i] * np.sqrt(Tini) # Assigning most probable velocities
                vy[i] = vy[i] * np.sqrt(Tini) # Assigning most probable velocities
                vz[i] = vz[i] * np.sqrt(Tini)
    #return x,y,vx,vy


def sign(L, dx):
    if (dx >= 0.0):
        return abs(L)  # if dx is positive, then dx = dx - L.
    else:  # dx < 0
        return -abs(L) # if dx is negative, then dx = dx + L.

def Forces(t, w, PE, PEorW):
    #invr2=0.0
    r2cut = 9.0
    PE = 0.0
    for i in range(0,Natm):
        fx[i][t] = 0 # to start with zero.
        fy[i][t] = 0 # to start with zero.
        fz[i][t] = 0  # to start with zero.
    for i in range(0,Natm-1):
        for j in range(i+1, Natm):
            dx = x[i] - x[j]
            dy = y[i] - y[j]
            dz = z[i] - z[j]
            if (abs(dx) > 0.5 * L):
                dx = dx - sign(L,dx)
            if (abs(dy) > 0.5 * L):
                dy = dy - sign(L,dy)
            if (abs(dz) > 0.5 * L):
                dz = dz - sign(L,dz)
            r2 = dx * dx + dy * dy + dz * dz

            if (r2 < r2cut):
                if (r2 == 0):
                    r2 = 0.0001
                invr2 = 1.0/r2

                wij = 48 * ((invr2**3) - 0.5) * (invr2**3) # Magnitude of Force
                #print(r2, invr2,wij)
                fijx = wij * dx * invr2
                fijy = wij * dy * invr2
                fijz = wij * dz * invr2
                fx[i][t] = fx[i][t] + fijx  # fx of ith particle
                fy[i][t] = fy[i][t] + fijy  # fy of ith particle
                fz[i][t] = fz[i][t] + fijz
                fx[j][t] = fx[j][t] - fijx  # fx of jth particle
                fy[j][t] = fy[j][t] - fijy  # fy of jth particle
                fz[j][t] = fz[j][t] - fijz
                # these fij = - fji are opposite of each other => Negative sign in one of them.


                PE = PE  + 4.0*(invr2**3)*((invr2**3) - 1.0)
                w = w + wij

    if (PEorW == 1):
        return PE
    else:
        return w


def Tevolution():
    f = open('KE.dat', 'w')
    g = open('PE.dat', 'w')
    hh = open('TE.dat', 'w')
    gg = open('temp.dat', 'w')
    avKE = 0
    avPE = 0
    avT = 0
    t1 = 0
    h = 0.031
    h_by2 = h/2.0
    KE = 0.0
    PE = 0.0
    w = 0.0
    initial_pos_vel()
    for i in range(0,Natm):
        KE = KE + (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i])*0.5 # KE calculation at the start up. These random velocities should
        # take up values following Maxwellian distribution.

    PE = Forces(t1, w, PE, 1) # forces are computed here for the first time when particles are placed on a 2D-grid.
    time = 1
    while time <= 1000:
        print(time)
        #print("hello")
        #t1 += 1
        for i in range(0,Natm):
            PE  = Forces(t1,w, PE, 1) # Calculate the forces here.
            x[i] = x[i] + h*(vx[i] + h_by2 * fx[i][t1]) # Velocities will be updated here in next iteration.
            if x[i] <= 0.0:  # PBC for x
                x[i] = x[i] + L
            if x[i] > L:
                x[i] = x[i] - L
            #print(x[i],'pos')
            y[i] = y[i] + h*(vy[i] + h_by2 * fy[i][t1])
            if y[i] <= 0.0: # PBC for y
                y[i] = y[i] + L
            if y[i] > L:
                y[i] = y[i] - L
            z[i] = z[i] + h * (vz[i] + h_by2 * fz[i][t1])
            if z[i] <= 0.0:  # PBC for y
                z[i] = z[i] + L
            if z[i] > L:
                z[i] = z[i] - L
        # Once this loop finishes, it generates a microstate and a set of forces on all particles of that microstate.
        t2 = 1
        PE = 0
        PE = Forces(t2, w , PE, 1) # Calculating new forces for previous microstate,
        # it is mandatory before computing velocities.
        KE = 0
        w = 0
        for i in range(0,Natm):
            vx[i] = vx[i] +  h_by2 * (fx[i][t1] + fx[i][t2]) # Updating velocities along x
            vy[i] = vy[i] + h_by2 * (fy[i][t1] + fy[i][t2]) # Updating velocities along y
            vz[i] = vz[i] + h_by2 * (fz[i][t1] + fz[i][t2])  # Updating velocities along z
            KE = KE + (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i])*0.5 # KE calculation.

    #if (t == 0):
    #    t = 1
        avKE = avKE + KE
        avPE = avPE + PE
        T = KE/Natm
        avT = avT + T
        time += 1
        t = time
        Tavg = avT / t
        avgKE = avKE / t
        avgPE = avPE / t

        f.write(str(time) + ' ' + str(avgKE) + '\n')
        g.write(str(time) + ' ' + str(avgPE) + '\n')
        hh.write(str(time) + ' ' + str(avgKE + avgPE) + '\n')
        gg.write(str(time) + ' ' + str(Tavg) + '\n')
Tevolution()
