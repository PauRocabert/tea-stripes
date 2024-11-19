from particle_class import _particle_
import numpy as np 
import matplotlib.pyplot as plt


N = 100 #particles
T = 373
m = 1
dt = 0.01
c = 0.01
time_steps = 10000

particles = [_particle_(T, m,c, n) for n in range(N)]
x = np.zeros(shape=(N,time_steps))
y = np.zeros(shape=(N,time_steps))
for t in range(time_steps):
    forces = [particle.force(particles) for particle in particles]
    for n,particle in enumerate(particles):
        particle.brownian_step(forces[n],dt)
        x[n][t], y[n][t] = particle.pos

np.savez('positions_brownian.npz', x=x, y=y)



        
  
    

        

        
    
        
