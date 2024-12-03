import numpy as np 

density = 2.12e3
m = 4/3*np.pi*(1e-4)**3*density 
mu = 8.9e-4
A0 = 1e-9
k = 2*np.pi/(5e-3)
g = 9.81
omega = 2*np.pi/0.1
a0 = m*g/(mu*omega)*k**2*A0
print(a0)