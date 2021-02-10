import numpy as np
import matplotlib.pyplot as plt

#x,y = np.loadtxt("xy_NACA64A010", delimiter=" ", unpack=True, skiprows=1, usecols=(2,3))
#x,y = np.loadtxt("xy_hayabusa", delimiter=" ", unpack=True, skiprows=1, usecols=(2,3))
#x,y = np.loadtxt("xy_NACA4412", delimiter=" ", unpack=True, skiprows=1, usecols=(2,3))

x,y = np.loadtxt("External shape", delimiter=" ", unpack=True, skiprows=1)

plt.plot(x,y,"o",color="blue", markersize=1)
#plt.xlim(-0.32,0.12)
#plt.ylim(-0.22,0.22)
#plt.xlim(-0.1,1.1)
#plt.ylim(-0.1,1.1)

plt.show()