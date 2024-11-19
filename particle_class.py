import numpy as np
import math



e = 1.6021766e-19
eta = 1e-3 #Pa
epsilon0 = 8.85541878e-12
epsilon_r = 80.10
epsilon = epsilon_r*epsilon0
Na = 6.022e23
kb = 1.3806488e-23
class _particle_():
    def __init__(self, T, m,c, index):
        I = c*Na*1e3
        kappa = 2.2912074e-3*np.sqrt(I/T)
        self.R = 1e-6 #m
        self.kappa = kappa
        surface_density = 6e-3 #estimar millor 
        self.DebyeHuckel = 2*np.pi*(surface_density**2)/(kappa*epsilon)
        self.VanderWaals = 1e-20*kappa**2/(12*np.pi) 
        sigma = np.sqrt(m/(2*math.pi*T))
        self.T = T
        self.D = kb*self.T/(6*np.pi*eta*self.R)
        self.m = m
        self.index = index
        self.v = np.array([np.random.normal(loc=0, scale = sigma), np.random.normal(loc=0, scale = sigma)])
        self.pos = np.array([np.random.uniform(low= -5000,high=5000),np.random.uniform(low=-5000, high=5000)])

    def force(self, particles):
        f = np.zeros(2)
        count = 0
        for particle in particles:
            if particle.index != self.index:
                dr = particle.pos - self.pos
                distance = np.sqrt(np.dot(dr,dr))
                if distance <self.R*self.kappa:
                    Colision = 10000
                    count +=1
                else:
                    Colision = 0
                if distance < 4*self.R*self.kappa:
                    VanderWaals_force = -self.VanderWaals/(distance**2)  
                else:
                    VanderWaals_force = 0
                Debye_Huckel_force = self.DebyeHuckel*np.exp(-distance)
                f+= (VanderWaals_force + Debye_Huckel_force + Colision)*dr/distance
        return f
    
    def brownian_step(self, force, dt):
        acc = (1/(np.sqrt(2*self.D))*1/(6*np.pi*eta))*force
        rdot = acc + np.random.normal(size=2)
        self.pos += rdot*dt
    def langevin_step(self, force, dt):
        zeta = 6*np.pi*self.R*eta
        d2rdt2 = 1/self.m*(force + np.sqrt(2*zeta*kb*self.T)*np.random.normal(size=2) - zeta*self.v)
        self.v += d2rdt2*dt
        self.pos = self.v*dt
    
    





        

