import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from euler_brownian import positions
from correlation import analysis
import os 

list_dir = os.listdir('./')
Temp = [300,320,320,340,360,380]
con = [0.001, 0.01, 0.1,1.0]
min_dist = np.zeros(shape=(len(Temp), len(con)))
for i,T in enumerate(Temp):
    for j,c in enumerate(con):
        data = f'positions_brownian_c_{c}_T_{T}.npz'
        if not(data in list_dir):
            print(data)
            positions(T=T, c=c, xy=True)
            min_dist[i][j]=analysis(data,c=c, t=T)
        if not(f'animation_c={c}_T={T}.mp4' in list_dir):
            x,y=np.load(data)['x'], np.load(data)['y']
            Na = 6.022e23
            I = c*Na*1e3
            radius = 1e-6
            kappa = 2.2912074e-3*np.sqrt(I/T)
            N, t = np.shape(x)
            fig, ax = plt.subplots()
            scater = ax.scatter(x[:,0], y[:,0], s=1)
            ax.set_xlim(np.mean(x)-3*np.std(x), np.mean(x)+3*np.std(x))
            ax.set_xlabel('x/r')
            ax.set_ylabel('y/r')
            ax.set_ylim(np.mean(y)-3*np.std(y), np.mean(y)+3*np.std(y))
            def update(frame):
                scater.set_offsets(np.c_[x[:, frame], y[:, frame]])
                return scater,
            ani = FuncAnimation(fig, update, frames=int(t), blit=True)
            ani.save(f'animation_c={c}_T={T}.mp4', fps=100) 
    fig, ax = plt.subplots()
    ax.plot(con,min_dist[i])
    ax.set_xlabel('c (M)')
    ax.set_ylabel('<d>/r')
    plt.savefig(f'min_dist_{T}.png', dpi=200)
    
np.savez(f'min_dist.npz', min_dist=min_dist)

        
