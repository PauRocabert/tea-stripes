import numpy as np
import math
from gauss_seidel import vel



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
        phi0 = 0.03
        self.R = 1e-4 #m
        self.kappa = kappa
#bulk DLVO

#modified DLVO
        z = np.tanh(e*phi0/(4*kb*T))
        self.DLVO = 64*c*kb*T*(z**2)/(self.R)
#Van der Waals
        self.VanderWaals = 0.95e-20*self.kappa**3/(12*np.pi)
        sigma = np.sqrt(m/(2*math.pi*T))
        self.T = T
#assumin partícules esfèriques
        self.D = kb*self.T/(6*np.pi*eta*self.R)
        self.m = m
        self.index = index
        self.lengh = 2*np.pi/(2.64)*1e-2*kappa #wavelenght 
#        self.v = np.array([np.random.normal(loc=0, scale = sigma), np.random.normal(loc=0, scale = sigma)])
#        self.pos = np.array([np.random.randint(low = 0, high = 2)*self.lengh+np.random.normal(loc=0, scale = self.lengh/20), np.random.normal(loc=0, scale = 10)])
        self.pos = np.array([ np.random.uniform(low= -self.R*kappa, high=self.R*kappa), np.random.uniform(low= -self.R*kappa, high=self.R*kappa)])*1000
        self.z_max = self.R*kappa*2000
        self.F0 = 20e5
        self.l = 1e-9

    def force_surface(self, particles):
        f = np.zeros(2)
        for particle in particles:
            if particle.index != self.index:
                dr = particle.pos - self.pos
                distance = np.sqrt(np.dot(dr,dr))
                VanderWaals_force = self.VanderWaals*(2 + 5.32/(100*self.kappa)*distance)/((distance**3)*(1+5.32/(100*self.kappa)*distance))               
                Debye_Huckel_force = self.DLVO*np.exp(-distance)
                hydrophobic_force = self.F0*np.exp(-distance/(self.kappa*self.l))
                f+= (-VanderWaals_force -hydrophobic_force+ Debye_Huckel_force)*dr/distance
        return f
    
    def force_vertical(self, particles):
        lambda_b = 1/(4*np.pi*kb*self.T*epsilon)*e**2
        Yukanawa_coef = lambda_b*self.kappa**2*(np.exp(self.kappa*self.R)/(1+self.kappa*self.R))**2
        f = np.zeros(2)
        for particle in particles:
            if particle.index != self.index:
                dr = particle.pos - self.pos
                distance = np.sqrt(np.dot(dr,dr))
                VanderWaals_force = self.VanderWaals*(2 + 5.32/(100*self.kappa)*distance)/((distance**3)*(1+5.32/(100*self.kappa)*distance))               
                Debye_Huckel_force = 2*Yukanawa_coef*np.exp(-distance)/(distance)
                hydrophobic_force = self.F0*np.exp(-distance/(self.kappa*self.l))
                buyancy = (1.0 - 2.12)*1000*9.81*(4/3*np.pi*(self.R)**3)*np.array([0,1]) #tannis acid density
                f+= (-VanderWaals_force -hydrophobic_force+ Debye_Huckel_force + buyancy)*dr/distance
        return f

    def brownian_step_xy(self, force, dt, convection):
        f = 1/(2*self.kappa*kb*self.T)*force
        vmax = 300e-5 #m/s
        print(f)
        rdot =  f + np.random.normal(scale = np.sqrt(dt), size=2) + convection*np.sin(self.pos[0]*2*np.pi/self.lengh)*(vmax/(2*self.D*self.kappa))
        self.pos += rdot*dt
    def brownian_step_xz(self, force, dt, convection):
        f = 1/(2*self.kappa*kb*self.T)*force
        rdot =  f + np.random.normal(scale = np.sqrt(dt), size=2) + convection*np.array([np.sin(self.pos[0]*2*np.pi/self.lengh),vel(self.pos[1]/self.z_max)])*(1/(2*self.D*self.kappa))
        self.pos += rdot*dt

    
    





        

