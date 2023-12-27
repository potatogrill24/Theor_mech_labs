import matplotlib.pyplot as p
from matplotlib.animation import FuncAnimation
import numpy as n

T = n.linspace(1, 2.5, 100)

Phi = n.sin(1.5*T) + 0.7

fgr = p.figure()
plt = fgr.add_subplot(1, 1, 1)
plt.axis('equal')

l = 1
r = 0.3

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
Sy = 1.5 - r - l * n.sin(Phi[0])

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
    Yb = 0.5 + r + l * n.cos(Phi[i])
    Sy = 1.5 - r - l * n.sin(Phi[i])
    Pruzhina.set_data(Xb + Xp, Sy + (Yb - Sy) * Yp)
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
