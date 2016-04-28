import numpy as np
import matplotlib.pyplot as plt
import warnings
from numba import jit

# ---------------------------------------------------- #
# Rotate shape around a point
@jit
def rotateShape(x, y, x0, y0, theta):
    theta = theta / 180 * np.pi
    x = x - x0
    y = y - y0
    x_ = np.zeros(x.shape)
    y_ = np.zeros(y.shape)
    for i in range(0, len(x)):
        x_[i] = np.cos(theta) * x[i] - np.sin(theta) * y[i]
        y_[i] = np.sin(theta) * x[i] + np.cos(theta) * y[i]
    # return (x, y)
    return (x_ + x0, y_ + y0)
# ---------------------------------------------------- #

# ---------------------------------------------------- #
# GENERATE DOBLE WEDGE SHAPE
@jit
def generateWedge(a, b, n, type, plotShape=False):
    # Define x and y coordinate
    xFrontTop_ = np.linspace(0, a / 2.0, n);                    xFrontTop = xFrontTop_
    xBackTop_ = np.linspace(a / 2.0, a, n);                     xBackTop = xBackTop_[1:]
    xBackBottom_ = np.flipud(np.linspace(a / 2.0, a, n));       xBackBottom = xBackBottom_[1:]
    xFrontBottom_ = np.flipud(np.linspace(0, a / 2.0, n));      xFrontBottom = xFrontBottom_[1:-1]

    yFrontTop_ = b / a * xFrontTop_;                            yFrontTop = yFrontTop_
    yBackTop_ = -b / a * (xBackTop_ - a / 2.0) + b / 2;         yBackTop = yBackTop_[1:]
    yBackBottom_ = b / a * (xBackBottom_ - a / 2.0) - b / 2;    yBackBottom = yBackBottom_[1:]
    yFrontBottom_ = -b / a * xFrontBottom_;                     yFrontBottom = yFrontBottom_[1:-1]

    if type == 'twoSided':
        # x = np.concatenate((xFrontTop, xBackTop, xBackBottom, xFrontBottom))
        # y = np.concatenate((yFrontTop, yBackTop, yBackBottom, yFrontBottom))
        xTop = np.concatenate((xFrontTop, xBackTop))
        xBottom = np.concatenate((xBackBottom, xFrontBottom))
        yTop = np.concatenate((yFrontTop, yBackTop))
        yBottom = np.concatenate((yBackBottom, yFrontBottom))
    elif type == 'oneSided':
        x = np.concatenate((xFrontTop, xBackTop))
        y = np.concatenate((yFrontTop, yBackTop))

    if plotShape:
        plt.figure()
        plt.plot(x, y, 'k')
        plt.xlim([-5, 14])
        plt.ylim([-5, 5])
        plt.show()

    return (x, y)
    # return (xTop, xBottom, yTop, yBottom)
# ---------------------------------------------------- #

# ---------------------------------------------------- #
# Generate Line
@jit
def generateLine(x0, y0, x1, y1, n):
    x = np.linspace(x0, x1, n)
    y = np.linspace(y0, y1, n)
    return (x, y)
# ---------------------------------------------------- #

# ---------------------------------------------------- #
# Calculate surface gradient for piston theory calculations
# Z(x, t) = Position of airfoil surface
# vn = Normal velocity of airfoil surface
#    = dZ/dt + V dZ/dx
@jit
def calcSurfaceGradient(x, y):
    # The following is added to suppress zero division warning
    warnings.simplefilter('ignore')

    x = np.concatenate((x, [x[0]]))
    y = np.concatenate((y, [y[0]]))
    dydx = (y[2:] - y[:-2]) / (x[2:] - x[:-2])
    dydx = np.concatenate(([(y[1] - y[-2]) / (x[1] - x[-2])], dydx))

    # Converts inf values to zero
    dydx[np.isinf(dydx)] = 0

    return dydx
# ---------------------------------------------------- #

# ---------------------------------------------------- #
# Calculate pressure at each node using piston theory
@jit
def calcSurfacePressure(x, y, ydot, pInfty=6.934, aInfty=343.0, gamma=1.4, V=1374.7716, rhoInfty=1.225):
    # gamma:    For diatomic gas
    # aInfty:   Speed of sound (m/s)
    # pInfty:   Free stream pressure
    # rhoInfty: Free stream density
    # V:        Free stream velocity
    dZdx = calcSurfaceGradient(x, y)
    dZdt = 0.0
    vn = dZdt + V * dZdx
    # P = pInfty * (1 + ((gamma - 1) / 2) * (vn / aInfty))**(2 * gamma / (gamma - 1))
    # P = pInfty + rhoInfty * aInfty**2 * ((vn / aInfty) +
    #                                      (gamma + 1) / 4 * (vn / aInfty)**2 +
    #                                      (gamma + 1) / 12 * (vn / aInfty)**3)
    P = pInfty + \
        pInfty * (gamma * vn / aInfty +
                  gamma * (gamma + 1) / 4 * (vn / aInfty)**2.0 +
                  gamma * (gamma + 1) / 12 * (vn / aInfty)**3.0)
    # P = np.divide(P, pInfty)
    return P
# ---------------------------------------------------- #

# ---------------------------------------------------- #
# Calculates tangent and normal vectors at each point on the surface
@jit
def calcSurfaceVectors(x, y):
    x = np.concatenate((x, [x[0]]))
    y = np.concatenate((y, [y[0]]))
    tangentVectorX = x[1:] - x[:-1]
    tangentVectorY = y[1:] - y[:-1]

    theta = np.pi/2
    rotationMatrix = np.matrix([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    normalVectorX = -tangentVectorY
    normalVectorY = tangentVectorX
    return (tangentVectorX, tangentVectorY, normalVectorX, normalVectorY)
# ---------------------------------------------------- #

# ---------------------------------------------------- #
# Calculates the forces in x and y directions at each point on the surface
@jit
def calcSurfaceForce(P, x, y):
    ds = np.sqrt((x[1] - x[0])**2.0 + (y[1] - y[0])**2.0)
    [tangentVectorX, tangentVectorY, normalVectorX, normalVectorY] = calcSurfaceVectors(x, y)
    Fx = -np.multiply(normalVectorX, P) * ds
    Fy = -np.multiply(normalVectorY, P) * ds
    return (Fx, Fy)
# ---------------------------------------------------- #