import matplotlib.pyplot as p
from matplotlib.animation import FuncAnimation
import numpy as n

from scipy.integrate import odeint

def SstDiffEq(y, t, m1, m2, m3, l, M, c, g):
    #y' = [phi, S, phi', S'] -> dy = [phi', S', phi'', S'']
    dy = n.zeros_like(y)
    dy[0] = y[2]
    dy[1] = y[3]

    #a11 * phi'' + a12 * S'' = b1
    #a21 * phi'' + a22 * S'' = b2

    a11 = 2 * ((3 * m1 * n.sin(y[0])**2) + (2 * m2 * n.cos(y[0])**2)) * l**2
    a12 = 0
    b1 = 4 * M - (3 * m1 - 2 * m2) * l**2 * y[2]**2 * n.sin(2*y[0]) - 4 * (c*(l * n.sin(y[0]) - y[1]) + m2 * g) * l * n.cos(y[0])

    a21 = 0
    a22 = m3
    b2 = c * (l * n.sin(y[0] - y[1])) - m3 * g

    dy[2] = (b1 * a22 - b2 * a12) / (a11 * a22 - a12 * a21)
    dy[3] = (b2 * a11 - b1* a21) / (a11 * a22 - a12 * a21)

    return dy

m1 = 2
m2 = 0.2
m3 = m2
l = 1
M = 3
c = 5
g = 9.81
r = 0.3

T = n.linspace(0, 10, 300)
y0 = [0, 0, n.pi/3, 0]

Y = odeint(SstDiffEq, y0, T, (m1, m2, m3, l, M, c, g))
print(Y)

Phi = Y[:,0] + 1.7
S = 1.5 - r - l * n.sin(Phi[0])
S = Y[:,1]
Phit = Y[:,2]
St = Y[:,3]

fgrp = p.figure()
plPhi = fgrp.add_subplot(4,1,1)
plPhi.plot(T, Phi)

plS = fgrp.add_subplot(4,1,2)
plS.plot(T, S)

Phitt = n.zeros_like(T)
Stt = n.zeros_like(T)

for i in range(len(T)):
    Phitt[i] = SstDiffEq(Y[i], T[i],  m1, m2, m3, l, M, c, g)[1]
    Stt[i] = SstDiffEq(Y[i], T[i],  m1, m2, m3, l, M, c, g)[2]

NX = -m1 * l * (Phitt * n.sin(Phi) + Phit**2 * n.cos(Phi)) / 2
NY = m2 * l * (Phitt * n.cos(Phi) + Phit**2 * n.sin(Phi)) + m3 * Stt + (m1 + m2 + m3) * g

plNX = fgrp.add_subplot(4,1,3)
plNX.plot(T, NX)

plNY = fgrp.add_subplot(4,1,4)
plNY.plot(T, NY)

# p.show()

# ===================anima=============

fgr = p.figure()
plt = fgr.add_subplot(1, 1, 1)
plt.axis('equal')

#оси
plt.plot([0, 0], [0, 3])
plt.plot([3, 0], [0, 0])
#стрелки для осей
plt.plot([0, -0.05], [3, 2.85])
plt.plot([0, 0.05], [3, 2.85])
plt.plot([3, 2.85], [0, 0.05])
plt.plot([3, 2.85], [0, -0.05])
#линии по бокам от оси y
plt.plot([-0.1, -0.1], [0, 2.7])
plt.plot([0.1, 0.1], [0, 2.7])

#шаблон окружности
Alp = n.linspace(0, 2*n.pi, 100)
Xc = r * n.cos(Alp)
Yc = r * n.sin(Alp)

Xa = l * n.sin(Phi[0])
Ya = 0.5 + r

Disk = plt.plot(Xc + Xa, Yc + Ya)[0]
Xb = 0
Yb = 1.5 + r + l * n.cos(Phi[0])
Sx = 0

AB = plt.plot([Xa, Xb], [Ya, Yb])[0]

#Шаблон пружины

Np = 30
Yp = n.linspace(0,1, 2*Np + 1)
Xp = 0.1 * n.sin(n.pi/2 * n.arange(2*Np + 1))

Pruzhina = plt.plot(Xb + Xp, Ya + (Yb - Ya) * Yp)[0]

#Шаблон спиральной пружины у диска

Ns = 3
r1 = 0.01
r2 = 0.1
numponts = n.linspace(0,1,30*Ns + 1)
Betas = numponts * (2*n.pi * Ns - Phi[0])
Xs = n.sin(Betas) * (r1 + (r2 - r1)*numponts)
Ys = n.cos(Betas) * (r1 + (r2 - r1)*numponts)

SpPruzhina = plt.plot(Xs + Xa, Ys + Ya)[0]

#Шаблон спиральной пружины у конца стержня

Ns = 3
r1_1 = 0.02
r2_1 = 0.1
numponts = n.linspace(0,1,30*Ns + 1)
Betas_1 = numponts * (2*n.pi * Ns - Phi[0])
Xs_1 = n.sin(Betas_1) * (r1 + (r2 - r1)*numponts)
Ys_1 = n.cos(Betas_1) * (r1 + (r2 - r1)*numponts)

SpPruzhina_1 = plt.plot(Xs_1 + Xb, Ys_1 + Yb)[0]

def run(i):
    Xa = l * n.sin(Phi[i])
    Disk.set_data(Xc + Xa, Yc + Ya)
    Yb = 1.5 + r + l * n.cos(Phi[i])
    S = 1.5 - r - l * n.sin(Phi[i])
    Pruzhina.set_data(Xb + Xp, S + (Yb - S) * Yp)
    AB.set_data([Xa, Xb], [Ya, Yb])

    
    Betas = numponts * (2*n.pi * Ns - Phi[i])
    Xs = n.sin(Betas) * (r1 + (r2 - r1)*numponts)
    Ys = n.cos(Betas) * (r1 + (r2 - r1)*numponts)
    SpPruzhina.set_data(Xs + Xa, Ys + Ya)
    
    Betas_1 = numponts * (2*n.pi * Ns - Phi[i])
    Xs_1 = n.sin(Betas_1) * (r1 + (r2 - r1)*numponts)
    Ys_1 = n.cos(Betas_1) * (r1 + (r2 - r1)*numponts)
    SpPruzhina_1.set_data(Xs_1 + Xb, Ys_1 + Yb)
    return

anim = FuncAnimation(fgr, run, frames = len(T), interval = 1)

p.show()
