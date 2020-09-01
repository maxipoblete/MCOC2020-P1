from matplotlib.pylab import *
from scipy.integrate import odeint
from numpy import *
import matplotlib.pyplot as plt
import warnings
import sys
warnings.simplefilter("ignore")


m = 1 
f = 1
ε = 0.2
w = 2*pi*f
k = m*w**2
c = 2*ε*w*m

β = c/(2*m)

z0 = array([1,1])

def solucion_odeint(z,t):
    zp = zeros(2)
    zp[0] = z[1]
    zp[1] = -(2*β*z[1]) - (z[0]*(w**2))
    return zp



def eulerint(solucion_odeint,z0,t,Nsubdivisiones):
    Nt = len(t)
    Ndim = len(array(z0))
    z = zeros((Nt,Ndim))
    
    z[0,:] = z0
    # Z_(i+1) = Zp_i * dt + z_i
    
    for i in range(1,Nt):
        t_anterior = t[i-1]
        dt = (t[i] - t[i-1])/float(Nsubdivisiones)
        z_temp = z[i-1,:].copy()
        for k in range(Nsubdivisiones):
            z_temp +=  dt * solucion_odeint(z_temp, t_anterior+k*dt)
        z[i,:] = z_temp
        
    return z
        

t = linspace(0,4,100)
z0 = array([1,1])
sol = odeint(solucion_odeint,z0,t)
z_odeint = sol[:,0]


z_euler1 = eulerint(solucion_odeint,z0,t,Nsubdivisiones=1)
z_euler10 = eulerint(solucion_odeint,z0,t,Nsubdivisiones=10)
z_euler100 = eulerint(solucion_odeint,z0,t,Nsubdivisiones=100)

plt.plot(t,z_euler1[:,0],"--",color="green",linewidth=1,label="Euler Nsubs = 1")
plt.plot(t,z_euler10[:,0],"--",color="red",linewidth=1,label="Euler Nsubs = 10")
plt.plot(t,z_euler100[:,0],"--",color="orange",linewidth=1,label="Euler Nsubs = 100")
plt.plot(t,z_odeint,color="blue",linewidth=1,label="Odeint")
plt.plot(t,((0.28461)**t)*cos(6.15624*t)+0.366561*((0.28461)**t)*sin(6.15624*t),color="k",linewidth=2,label="Analitica")
plt.legend()
plt.xlabel("Tiempo t (segundos)")
plt.ylabel("Distancia x(t) (metros)")
plt.show()






