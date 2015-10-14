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

gamma = 1.4 # For diatomic gas
ainfty = 343.0 # Speed of sound (m/s)
Pinfty = 10.0 # Free stream pressure
V = 10.0 # Free stream velocity
# airfoil = np.loadtxt('NACA0012.dat')
airfoil = np.loadtxt('NACA2414.dat')
# airfoil = np.loadtxt('NACA63-412.dat')

# Nomenclature
# Z(x, t) = Position of airfoil surface
# vn = Normal velocity of airfoil surface
#    = dZ/dt + V dZ/dx

# Define the surface of airfoil as part of a circle
P1 = [-1, 0]
P2 = [1, 0]
R = 1.5
def equation_of_circle(X, p1=P1, p2=P2, r=R):
    xc, yc = X
    x1, y1 = p1
    x2, y2 = p2
    return ((x1 - xc)**2 + (y1 - yc)**2 - r**2, (x2 - xc)**2 + (y2 - yc)**2 - r**2)
xc, yc = fsolve(equation_of_circle, (0, -1))
theta = np.arcsin(np.abs(P2[0] - P1[0]) / 2 / R)
theta = np.linspace(np.pi / 2 + theta, np.pi / 2 - theta, 100)
Zx = xc + R * np.cos(theta)
Zy = yc + R * np.sin(theta)

# Comment if you don't want to use airfoil
Zx = airfoil[:, 0]
Zy = airfoil[:, 1]

# plt.figure()
# plt.plot(Zx, Zy)
# plt.xlim([P1[0], P2[0]])
# plt.axis('equal')
# plt.show()

# Calcualte the spatial part of the velocity
dZdx = (Zy[2:] - Zy[:-2]) / (Zx[2:] - Zx[:-2])
# Calcualte the temporal part of the velocity
dZdt = 0.0
# Calculate normal velocity
vn = V * dZdx

# Calculate pressure
P = Pinfty + Pinfty * (gamma * vn / ainfty + gamma * (gamma + 1) * (vn / ainfty)**2 / 4 +
                       gamma * (gamma + 1) * (vn / ainfty)**3 / 12)

# plt.figure()
# plt.plot(Zx[1:-1], P)
# plt.xlabel('Location on X axis (m)')
# plt.ylabel('Pressure (Pa)')
# plt.show()

fig, ax1 = plt.subplots()
ax1.plot(Zx, Zy, 'k.')
ax1.set_xlabel('Location on X axis (m)')
# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel('Location on Y axis (m)', color='k')
plt.axis('equal')
for tl in ax1.get_yticklabels():
    tl.set_color('k')

ax2 = ax1.twinx()
ax2.plot(Zx[1:-1], P, 'ro')
ax2.set_ylabel('Pressure', color='r')
for tl in ax2.get_yticklabels():
    tl.set_color('r')
plt.show()
