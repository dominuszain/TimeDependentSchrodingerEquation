import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

data = np.load("wavepacket.npz")
frames = len(data["x1"])

fig, ax = plt.subplots(1, 1, figsize=(12, 12))
ax.set_facecolor("k")
ax.get_xaxis().set_ticks([])
ax.get_yaxis().set_ticks([])
ax.set_ylim(-2, 2)
ax.set_xlim(-1, 11)

points = [ax.plot(0, 0, marker="o")[0] for _ in range(11)]


def animate(i):
    for j, point in enumerate(points, start=1):
        point.set_data(data[f"x{j}"][i], data[f"y{j}"][i])
    return points


ani = animation.FuncAnimation(fig, animate, frames=frames, interval=50, blit=True)
ani.save("output.gif", writer="pillow", fps=50)
