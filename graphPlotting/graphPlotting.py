"""
Code to plot out animated 1D data from text file

@author: Gabriel Canet Espinosa
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

data = np.loadtxt('densities.txt', delimiter=',')
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)

def animate(i):
    x = data[i, ::2]
    y = data[i,1::2]
    ax1.clear()
    ax1.plot(x, y, 'o', mfc='none', ms=2)

ani = animation.FuncAnimation(fig, animate, frames=100, interval=50)

plt.show()