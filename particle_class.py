import numpy as np
import math


R = 1e-5 #m
e = 1.6021766e-19
eta = 1e-3 #Pa
epsilon0 = 8.85541878e-12
epsilon_r = 80.10
epsilon = epsilon_r*epsilon0
H_131 = 1 #Hammacker connstant
VdWaals = H_131*R/12
Na = 6.022e24
kb = 1.3806488e-23
class _particle_():
    def __init__(self, T, m, L,c, index):
        I = c*Na
        kappa = 2.2912074e-3*np.sqrt(I/T)
        self.kappa = kappa
        surface_density = 6e-3 #estudiar millor
        self.DebyeHuckel = 2*np.pi*(surface_density**2)/(kappa*epsilon)
        self.VanderWaals = 1e-21*kappa**2/(12*np.pi) 
        sigma = np.sqrt(m/(2*math.pi*T))
        self.T = T
        self.D = kb*self.T/(6*np.pi*eta*R)
        self.m = m
        self.index = index
        self.v = np.array([np.random.normal(loc=0, scale = sigma), np.random.normal(loc=0, scale = sigma)])
        self.pos = np.array([np.random.normal(loc = 0, scale = 0.1),np.random.normal(loc=0, scale = 0.1)])
    def force(self, particles):
        f = np.zeros(2)
        for particle in particles:
            if particle.index != self.index:
                dr = particle.pos - self.pos
                distance = np.sqrt(np.dot(dr,dr))
                if distance <R:
                    Colision = 10000
                else:
                    Colision = 0
                if distance < 3*R:
                    VanderWaals_force = -VdWaals/(distance**2)  
                else:
                    VanderWaals_force = 0
                Debye_Huckel_force = self.DebyeHuckel*self.kappa*np.exp(-self.kappa*distance)
                f+= (VanderWaals_force + Debye_Huckel_force+ Colision)*dr/distance
        return f
    
    def brownian_step(self, force, dt):
        rdot = (1/(6*np.pi*eta))*force + np.sqrt(2*self.D)*np.random.normal(size=2)
        self.pos += rdot*dt
    def langevin_step(self, force, dt):
        zeta = 6*np.pi*R*eta
        d2rdt2 = 1/self.m*(force) + np.sqrt(2*zeta*kb*self.T)*np.random.normal(size=2) - zeta*self.v
        self.v += d2rdt2*dt
        self.pos = self.v*dt





        

