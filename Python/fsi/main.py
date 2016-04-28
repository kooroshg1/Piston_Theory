__author__ = 'koorosh'
'''
Piston theory is an inviscid unsteady aerodynamic method used extensively in hypersonic aeroelasticity,
which predicts a point function relationship between the local pressure on a lifting surface and the normal component
of fluid velocity produced by the lifting surface motion. Here we use the third order expansion of "simple wave"
expression for the pressure on a piston.

This code is based on the following paper:
@article{thuruthimattam2002aeroelasticity,
  title={Aeroelasticity of a generic hypersonic vehicle},
  author={Thuruthimattam, BJ and Friedmann, PP and McNamara, JJ and Powell, KG},
  journal={AIAA Paper},
  volume={1209},
  pages={2002},
  year={2002}
}
'''

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import pdb
import pistonSolver
import time

gamma = 1.4 # For diatomic gas
aInfty = 343.0 # Speed of sound (m/s)
pInfty = 6.934 # Free stream pressure
rhoInfty = 1.225 # Free stream density
V = 1374.7716 # Free stream velocity

t = time.time()
[x, y] = pistonSolver.generateWedge(2.54, 0.254, 5, 'oneSided', plotShape=False)
[x, y] = pistonSolver.rotateShape(x, y, np.mean(x), np.mean(y), 2)
P = pistonSolver.calcSurfacePressure(x, y, ydot=0)
[Fx, Fy] = pistonSolver.calcSurfaceForce(P, x, y)
elapsed = time.time() - t
print(elapsed)

writeData = np.zeros([len(x), 4])
writeData[:, 0] = x
writeData[:, 1] = P
writeData[:, 2] = Fx
writeData[:, 3] = Fy
np.savetxt('pressure.txt', writeData,
           fmt='%-14.4f %-14.4e %-14.4e %-14.4e',
           delimiter=' ',
           header='x            Pressure       Fx             Fy')

print(sum(Fx))
print(sum(Fy))

[tangentVectorX, tangentVectorY, normalVectorX, normalVectorY] = pistonSolver.calcSurfaceVectors(x, y)

# plt.figure()
# plt.quiver(x, y, normalVectorX, normalVectorY)
# # plt.quiver(x, y, tangentVectorX, tangentVectorY)
# plt.hold('on')
# plt.plot(x, y, 'ko-', ms=10)
# plt.axis('equal')
# plt.show()

# plt.figure()
# plt.plot(x, Fx, 'ko', ms=10)
# plt.xlabel('x')
# plt.ylabel('$F_x$')
# # plt.ylim([-200, 200])
# plt.figure()
# plt.plot(x, y, 'ko-', ms=10)
# plt.axis('equal')
# plt.xlabel('x')
# plt.ylabel('Pressure')
# plt.show()

# Cp = (P - Pinfty) / (0.5 * V**2.0)
# plt.figure()
# plt.plot(Cp)
# plt.show()

# plt.figure()
# plt.plot(x[:-1], P, 'o')
# plt.xlabel('Location on X axis (m)')
# plt.ylabel('Pressure (Pa)')
# plt.show()

# fig, ax1 = plt.subplots()
# ax1.plot(Zx, Zy, 'k.')
# ax1.set_xlabel('Location on X axis (m)')
# # Make the y-axis label and tick labels match the line color.
# ax1.set_ylabel('Location on Y axis (m)', color='k')
# plt.axis('equal')
# for tl in ax1.get_yticklabels():
#     tl.set_color('k')
#
# ax2 = ax1.twinx()
# ax2.plot(Zx[1:-1], P, 'ro')
# ax2.set_ylabel('Pressure', color='r')
# for tl in ax2.get_yticklabels():
#     tl.set_color('r')
# plt.show()


