import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

x1 = np.loadtxt('x1.csv')
y1 = np.loadtxt('y1.csv')
x2 = np.loadtxt('x2.csv')
y2 = np.loadtxt('y2.csv')
x3 = np.loadtxt('x3.csv')
y3 = np.loadtxt('y3.csv')
x4 = np.loadtxt('x4.csv')
y4 = np.loadtxt('y4.csv')
x5 = np.loadtxt('x5.csv')
y5 = np.loadtxt('y5.csv')
x6 = np.loadtxt('x6.csv')
y6 = np.loadtxt('y6.csv')
x7 = np.loadtxt('x7.csv')
y7 = np.loadtxt('y7.csv')
x8 = np.loadtxt('x8.csv')
y8 = np.loadtxt('y8.csv')
x9 = np.loadtxt('x9.csv')
y9 = np.loadtxt('y9.csv')
x10 = np.loadtxt('x10.csv')
y10 = np.loadtxt('y10.csv')
x11 = np.loadtxt('x11.csv')
y11 = np.loadtxt('y11.csv')

def animate(i):
    #ln1.set_data([0, x1[i], x2[i]], [0, y1[i], y2[i]])
    point1.set_data(x1[i],y1[i])
    point2.set_data(x2[i],y2[i])
    point3.set_data(x3[i],y3[i])
    point4.set_data(x4[i],y4[i])
    point5.set_data(x5[i],y5[i])
    point6.set_data(x6[i],y6[i])
    point7.set_data(x7[i],y7[i])
    point8.set_data(x8[i],y8[i])
    point9.set_data(x9[i],y9[i])
    point10.set_data(x10[i],y10[i])
    point11.set_data(x11[i],y11[i])
    
fig, ax = plt.subplots(1,1, figsize=(12,12))
ax.set_facecolor('k')
ax.get_xaxis().set_ticks([])    # enable this to hide x axis ticks
ax.get_yaxis().set_ticks([])    # enable this to hide y axis ticks
#ln1, = plt.plot([], [], 'ro--', lw=3, markersize=8)
point1, = ax.plot(0,0, marker="o")
point2, = ax.plot(0,0, marker="o")
point3, = ax.plot(0,0, marker="o")
point4, = ax.plot(0,0, marker="o")
point5, = ax.plot(0,0, marker="o")
point6, = ax.plot(0,0, marker="o")
point7, = ax.plot(0,0, marker="o")
point8, = ax.plot(0,0, marker="o")
point9, = ax.plot(0,0, marker="o")
point10, = ax.plot(0,0, marker="o")
point11, = ax.plot(0,0, marker="o")
ax.set_ylim(-2,2)
ax.set_xlim(-1,11)
ani = animation.FuncAnimation(fig, animate, frames=6000, interval=50)
ani.save('output.gif',writer='pillow',fps=50)