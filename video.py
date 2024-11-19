import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

x,y=np.load('positions.npz')['x'], np.load('positions.npz')['y']
N, T = np.shape(x)

fig, ax = plt.subplots()

scater = ax.scatter(x[:,0], y[:,0], s=1)
ax.set_xlim(np.min(x), np.max(x))
ax.set_ylim(np.min(y), np.max(y))
def update(frame):
    scater.set_offsets(np.c_[x[:, frame], y[:, frame]])
    return scater,
ani = FuncAnimation(fig, update, frames=T, interval=100, blit=True)
ani.save('animation_1.mp4', fps=100)  # Uncomment to save as an mp4 file
plt.show()
