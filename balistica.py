import scipy as sp
from scipy.integrate import odeint
import warnings
warnings.simplefilter("ignore")

"""
===================================
       DEFINO LAS UNIDADES
===================================
"""
cm = 0.01                       # LA UNIDAD BASE ES EL METRO
inch = 2.54*cm                  # PARA TRANSFORMAR LAS UNIDADES DE LA PELOTA DE BOWLING
g = 9.81                        # LA ACELERACION DE GRAVEDAD EN LA TIERRA


"""
===================================
  DEFINO BOLA Y EL COEF ARRASTRE
===================================
"""

p = 1.225                       # DENSIDAD DEL AIRE EN KG/M**3
cd = 0.47                       # DRAG COEF PARA CUERPO ESFERICO
D = 8.5*inch                    # DIAMETRO DE PELOTA DE BOWLING
r = D/2                         # RADIO DE PELOTA DE BOWLING
A = sp.pi*r**2                  # SUPERFICIE CIRCULAR DE CARA DE LA BOLA
CD = 0.5*p*cd*A                 # COEFICIENTE DE ARRASTRE (GRANDE) 
m = 15                          # MASA DE LA BOLA DE PROYECTIL EN KG


"""
===================================
       DEFINO LA FUNCION
===================================

z : es el vector de estado 
z = (x,y,vx,vy)

"""

V = 0

def bala(z,t):
    zp = sp.zeros(4)               # SE CREA EL VECTOR "z PUNTO"
    zp[0] = z[2]                   # LA DERIVADA DE LA POSICION EN X ES LA VELOCIDAD VX
    zp[1] = z[3]                   # LA DERIVADA DE LA POSICION EN Y ES LA VELOCIDAD VY
    
    v = z[2:4]                     # SACAMOS LOS ULTIMOS DOS COMPONENTES DE Z, LAS VELOCIDADES v = (vx,vy)
    v[0] = v[0] - V                # APLICAMOS LA OPOSICION DEL VIENTO
    v2=sp.dot(v,v)                 # 
    vnorm = sp.sqrt(v2)            # ESTE Y EL PASO ANTERIOR ES PARA NORMALIZAR LA VELOCIDAD
    
    FD = -CD*v2*(v/vnorm)          # LA DRAG FORCE SERA IGUAL A MENOS UN VECTOR IGUAL A LA VELOCIDAD PERO 
                                   # MULTIPLICADO POR UN PARAMETRO PROPOCIONAL CON EL COEFICIENTE DE ARRASTRE 
                                   # Y EL SIGNO NEGATIVO ES POR QUE ES OPUESTO
    
    zp[2] = FD[0]/m                # LA "ACELERACION EN X"
    zp[3] = FD[1]/m - g            # LA "ACELERACION EN Y"
     
    return zp                      # RETORNAMOS z PUNTO DEL ESTADO O ESPACIO DE ESTADO





"""
===================================
       HAGO EL GRAFICO 
       SEGUN LO PEDIDO
===================================

El proyectil se lanza desde el origen
y con igual velocidad en x e y.
 
"""

t = sp.linspace(0,30,1001)          # DEFINIMOS EL TIEMPO
vi= 100*(1000./3600)                # DEFINIMOS LA VELOCIDAD vx0=vy0 INICIAL
z0 = sp.array([0,0,vi,vi])          # DEFINIMOS EL ESTADO INICIAL

import matplotlib.pyplot as plt
plt.figure(1,figsize=(8,8))
plt.xlim(0,150)
plt.ylim(0,50)
plt.xlabel(" X (m) ",fontsize=15,fontweight="bold")
plt.ylabel(" Y (m) ",fontsize=15,fontweight="bold")
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.title(" Trayectoria para Distintos Vientos \n",fontsize=20,fontweight="bold")
plt.grid(True)
linewidth = 4
alpha = 0.8


V=0
sol1 = odeint(bala,z0,t)  
x1 = sol1[:,0]
y1 = sol1[:,1]
plt.plot(x1,y1,"skyblue",linewidth=linewidth,alpha=alpha,label="V = 0 m/s")

V=10
sol2 = odeint(bala,z0,t) 
x2 = sol2[:,0]
y2 = sol2[:,1]
plt.plot(x2,y2,"orange",linewidth=linewidth,alpha=alpha,label="V = 10.0 m/s")

V=20
sol3 = odeint(bala,z0,t) 
x3 = sol3[:,0]
y3 = sol3[:,1]
plt.plot(x3,y3,"greenyellow",linewidth=linewidth,alpha=alpha,label="V = 20.0 m/s")
plt.legend(fancybox=True,shadow=True,fontsize="xx-large")

plt.hlines(0.2,0,150,"yellowgreen",linewidth=4) # ESTE ES EL PASTO
plt.hlines(0.1,0,150,"saddlebrown",linewidth=3) # ESTA ES LA TIERRA
plt.tight_layout()
plt.savefig('Proyectiles.png', dpi=500)
plt.close()












