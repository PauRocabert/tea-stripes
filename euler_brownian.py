from particle_class import _particle_
import numpy as np 
import matplotlib.pyplot as plt
import time 

N = 100 #particles
T = 300
m = 1
L = 0.01
dt = 0.001
c = 1
time_steps = 10000

particles = [_particle_(T, m, L,c, n) for n in range(N)]
breakpoint()
x = np.zeros(shape=(N,time_steps))
y = np.zeros(shape=(N,time_steps))
for t in range(time_steps):
    forces = [particle.force(particles) for particle in particles]
    for n,particle in enumerate(particles):
        particle.brownian_step(forces[n],dt)
        x[n][t], y[n][t] = particle.pos
np.savez('positions.npz', x=x, y=y)



        
  
    

        

        
    
        
