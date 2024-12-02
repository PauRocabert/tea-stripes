import numpy as np 
from scipy.stats import linregress
import matplotlib.pyplot as plt
def analysis(data, c, t):
    x,y=np.load(data)['x'][:,-1], np.load(data)['y'][:,-1]
    N = np.shape(x)[0]
    Na = 6.022e23
    I = c*Na*1e3
    kappa = 2.2912074e-3*np.sqrt(I/t)
    distance_matrix = np.zeros(shape=(N,N))

    for n in range(N):
        for m in range(N):
            if n<m:
                distance_matrix[n][m] = np.sqrt((x[m]-x[n])**2+(y[m]-y[n])**2)
                distance_matrix[m][n] =  distance_matrix[n][m]
    dr = 10
    rs = np.arange(0,int(np.max(distance_matrix)),dr)
    Rs = np.array([r*np.ones(100) for r in rs])

    m = 0
    min_dist = np.average(np.array([np.mean(np.sort(row)[:3]) for row in distance_matrix]))
    print(min_dist)

    for n in range(N):
        ac = np.zeros(shape=len(rs))
        for i in range(len(rs)-1):
            nums =np.sum((distance_matrix[n] > Rs[i])*(distance_matrix[n] < Rs[i+1]))
            ac[i] += nums/(2*rs[i]*dr + dr**2)
        m += linregress(rs, ac).slope
    ac,rs=ac[ac>0], rs[ac>0]
    print(ac)
    fig, ax = plt.subplots(figsize=(12,8))
    ax.plot(rs, ac)
    plt.savefig(f'g(r,c={c},T={t}).png', dpi=200)
    return min_dist