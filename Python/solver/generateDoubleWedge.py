__author__ = 'koorosh'
import numpy as np
import matplotlib.pyplot as plt

a = 2.54;
b = 0.254;
n = 40;

xFrontTop = np.linspace(0, a / 2.0, n).reshape(-1, 1)
xBackTop = np.linspace(a / 2.0, a, n).reshape(-1, 1)
# xBackBottom = np.linspace(a / 2.0, a, n).reshape(-1, 1)
# xFrontBottom = np.linspace(0, a / 2.0, n).reshape(-1, 1)
xBackBottom = np.flipud(np.linspace(a / 2.0, a, n).reshape(-1, 1))
xFrontBottom = np.flipud(np.linspace(0, a / 2.0, n).reshape(-1, 1))

yFrontTop = b / a * xFrontTop
yBackTop = -b / a * (xBackTop[1:] - a / 2.0) + b / 2
yBackBottom = b / a * (xBackBottom[1:] - a / 2.0) - b / 2
yFrontBottom = -b / a * xFrontBottom[1:-1]

# Generate only top surface
wedgeCoordinate = np.zeros([yFrontTop.shape[0] + yBackTop.shape[0], 2])
# Generate full shape
# wedgeCoordinate = np.zeros([yFrontTop.shape[0] + yBackTop.shape[0] + yBackBottom.shape[0] + yFrontBottom.shape[0], 2])

wedgeCoordinate[0:yFrontTop.shape[0], 0] = xFrontTop[:, 0]
wedgeCoordinate[yFrontTop.shape[0]:yFrontTop.shape[0] + yBackTop.shape[0], 0] = xBackTop[1:, 0]
# wedgeCoordinate[yFrontTop.shape[0] + yBackTop.shape[0]:yFrontTop.shape[0] + yBackTop.shape[0] + yBackBottom.shape[0], 0] = xBackBottom[1:, 0]
# wedgeCoordinate[yFrontTop.shape[0] + yBackTop.shape[0] + yBackBottom.shape[0]:, 0] = xFrontBottom[1:-1, 0]

wedgeCoordinate[0:yFrontTop.shape[0], 1] = yFrontTop[:, 0]
wedgeCoordinate[yFrontTop.shape[0]:yFrontTop.shape[0] + yBackTop.shape[0], 1] = yBackTop[:, 0]
# wedgeCoordinate[yFrontTop.shape[0] + yBackTop.shape[0]:yFrontTop.shape[0] + yBackTop.shape[0] + yBackBottom.shape[0], 1] = yBackBottom[:, 0]
# wedgeCoordinate[yFrontTop.shape[0] + yBackTop.shape[0] + yBackBottom.shape[0]:, 1] = yFrontBottom[:, 0]
#
# wedgeCoordinate[0:n, 1] = yFrontTop[:, 0]
# wedgeCoordinate[n:2 * n, 1] = yBackTop[:, 0]
# wedgeCoordinate[2 * n:3 * n, 1] = yBackBottom[:, 0]
# wedgeCoordinate[3 * n:4 * n, 1] = yFrontBottom[:, 0]

# print(wedgeCoordinate)

plt.figure()
plt.plot(wedgeCoordinate[:, 0], wedgeCoordinate[:, 1], 'k')
plt.xlim([-5, 14])
plt.ylim([-5, 5])
plt.show()

np.savetxt('double_wedge.txt', wedgeCoordinate)