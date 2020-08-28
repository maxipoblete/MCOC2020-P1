import scipy as sp
import numpy as np
from scipy.integrate import odeint
import warnings
warnings.simplefilter("ignore")


"""
===============================================================
            DEFINO LOS PARAMETROS INICIALES
===============================================================
"""

R = 6371000             # RADIO DE LA TIERRA                – [ m ]
r = R + (700000)        # DISTANCIA TIERRA A SATELITE       – [ m ]
G = 6.674e-11           # CONSTANTE GRAVITACIONAL           – [ N m**2 / kg**2]
M = 5.972e24            # MASA DE LA TIERRA                 – [ kg ]
w = 7.272205e-5         # VELOCIDAD ANGULAR DE LA TIERRA    – [ rad/s ]
atmosfera = 80000       # ALTURA DE ATMOSFERA EN            – [ m ]
v0_x = 0 / 3.6          # VELOCIDAD INICIAL EN X            – [ m/s ] 
v0_y = 6335.5           # VELOCIDAD INICIAL EN Y            – [ m/s ] 
v0_z = 0                # VELOCIDAD INICIAL EN Z            – [ m/s ] 



"""
===============================================================
         ESTABLEZCO LA FUNCION PARA LAS ORBITAS
===============================================================
"""

def satelite(z,t):
    
    
    # 1. DEFINO LOS PARAMETROS Y MATRICES DE ROTACION A UTILIZAR

    θ = (w*t)
    c = ((G*M)/(r**3))
    
    R   = sp.array([[    np.cos(θ)      ,   -np.sin(θ)    ,0     ],
                    [    np.sin(θ)      ,    np.cos(θ)    ,0         ],
                    [0                  ,0                ,1]     ])
    
    
    Rp2  = 2 * sp.array([[    -np.sin(θ)      ,    -np.cos(θ)  ,0   ],
                         [     np.cos(θ)      ,    -np.sin(θ)  ,0     ],
                         [0                  ,0               ,0]     ])
    
    
    Rpp = sp.array([      [    -np.cos(θ)      ,    np.sin(θ)     ,0    ] 
                    ,     [    -np.sin(θ)      ,   -np.cos(θ)    ,0     ],
                          [0                  ,0                 ,0]     ])
    





    # 2. DEFINO LOS ELEMENTOS DEL VECTOR DE ESTADO Z
    z1 = sp.array([z[0],z[1],z[2] ])
    z2 = sp.array([z[3],z[4],z[5]])
    # 3. AQUI CREO LA TRANSFORMACION PARA EL VECTOR DE ESTADO ZP = [ZP1 , ZP2]
    zp = sp.zeros(6)
    
    # --- ZP1
    zp[0] = z[3] 
    zp[1] = z[4] 
    zp[2] = z[5] 
    
    
    # --- ZP2
    z2p = -(c*z1) - (   (R.T)@((w**2)*Rpp)@(z1) ) - (   (R.T)@(w*Rp2)@(z2) )
    zp[3] = z2p[0] 
    zp[4] = z2p[1] 
    zp[5] = z2p[2] 
    return zp



"""
===============================================================
            IMPORTO LAS LIBRERIAS PARA GRAFICAR
===============================================================
"""



import matplotlib.pyplot as plt


plt.figure(1,figsize=(8,10))

"----  GRAFICO LA ORBITA  ----"
t = sp.linspace(0,12750,1001)
z0  = sp.array([r,0,0,v0_x,v0_y,v0_z])

sol = odeint(satelite,z0,t)
x = sol[:,0]
y = sol[:,1]
z = sol[:,2]

plt.plot(x,y,"orange",linewidth=1.5,label=f"V = {v0_y} m/s ")
plt.plot(r,0,"bo",markersize=5,label="Inicio")

"----  GRAFICO LA TIERRA  ----"
α = sp.linspace(0,360,360)
plt.plot(R*np.cos(np.radians(α)),R*np.sin(np.radians(α)),"c")


Rplot = R
for i in range(50):
    Rplot=Rplot-5700*i
    plt.plot(Rplot*np.cos(np.radians(α)),Rplot*np.sin(np.radians(α)),"c",alpha=0.2)




"----  GRAFICO LA ATMOSFERA  ----"
plt.plot((R+atmosfera)*np.cos(np.radians(α)),(R+atmosfera)*np.sin(np.radians(α)),"--r",linewidth=1,label="Atmosfera")



"""----  GRAFICO LOS EJES PARA LA PROPORCIONALIDAD 
                Y LIMITE DE RADIO DE TIERRA          ----"""

plt.vlines(0,-R,R)
plt.hlines(0,-R,R)





"----  PARAMETROS EXTRAS PARA GRAFICAR   ----"
hfont = {'fontname':'Avenir'}
plt.title("MCOCP1-E2 MAXIMILIANO POBLETE\n ORBITA SATELITAL",**hfont,fontsize=20)
plt.xticks(**hfont,rotation=90,fontsize=15)
plt.xticks(**hfont,fontsize=16)
plt.yticks(**hfont,fontsize=16)
plt.xlabel("\nPosicion en X (m)",**hfont,fontsize=20)
plt.ylabel("\nPosicion en Y (m)",**hfont,fontsize=20)
plt.ticklabel_format(useOffset=False, style='plain')
plt.grid(True)
plt.tight_layout()
plt.legend()
# plt.savefig("Grafico 0 P1.png",dpi=500)
plt.show()


"---- escritura de archivo   ----"
# archivo = open("salida.txt","w")
# for i in range(len(x)):
#     archivo.write(f"{x[i]} {y[i]} {z[i]}\n")
# archivo.close()
# print (x)
















#=========================================================================

#               EXTRAS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

#=========================================================================





#------------------------------------------------------------------------------
# Esto use para encontrar la velocidad limite
#------------------------------------------------------------------------------

# solucion = sp.linspace(6335.27,6335.29,100)
# radioactual = []

# t = sp.linspace(0,172800,1001)
# z0  = sp.array([r,0,0,v0_x,v0_y,v0_z])

# sol = odeint(satelite,z0,t)
# x = sol[:,0]
# y = sol[:,1]
# z = sol[:,2]

# for i in range(len(x)):
#     radioactual.append(np.sqrt (x[i]**2 + y[i]**2))
# min(radioactual)
# radiolimite = 6451000
# for s in solucion:
#     t = sp.linspace(0,14400,1001)
#     z0  = sp.array([r,0,0,v0_x,s,v0_z])
#     sol = odeint(satelite,z0,t)
#     x = sol[:,0]
#     y = sol[:,1]
#     z = sol[:,2]
#     radioactual = []
#     for i in range(len(x)):
#         radioactual.append(   np.sqrt (   x[i]**2    + y[i]**2  )    ) 
#     if radiolimite < min(radioactual):
#         print (f"todo ok {s}")     
#     if radiolimite >= min(radioactual):
#         print (min(radioactual))
#         print (f"problema {s}")
        

#------------------------------------------------------------------------------
# Esto use para graficar tiempo versus radio orbital
#------------------------------------------------------------------------------

# import matplotlib.pyplot as plt
# plt.figure(1,figsize=(8,10))
# t = sp.linspace(0,12750,1001)
# z0  = sp.array([r,0,0,v0_x,v0_y,v0_z])
# sol = odeint(satelite,z0,t)
# x = sol[:,0]
# y = sol[:,1]
# z = sol[:,2]
# radioactual=[]
# for i in range(len(x)):
#     radioactual.append(np.sqrt (x[i]**2 + y[i]**2))
# plt.plot(t/3600,np.array(radioactual),"orange",linewidth=1.5,label=f" Tiempo v/s Radio Orbital [V={v0_y} m/s] ")
# plt.ticklabel_format(useOffset=False, style='plain')
# plt.grid(True)
# plt.hlines(R,0,max(t)/3600,color="limegreen",linewidth=4,label="Superficie Terrestre - 6371 km")
# plt.hlines(R+atmosfera,0,max(t)/3600,"r","--",linewidth=1,label="Atmosfera (MESOPAUSA) - 6451 km")
# plt.hlines(6200000,0,max(t)/3600,color="peru",linewidth=100,alpha=0.2)
# hfont = {'fontname':'Avenir'}
# plt.title("MCOCP1-E2 MAXIMILIANO POBLETE\n TIEMPO V/S RADIO ORBITAL",**hfont,fontsize=20)
# plt.xticks(**hfont,rotation=90,fontsize=15)
# plt.xticks(**hfont,fontsize=16)
# plt.yticks(**hfont,fontsize=16)
# plt.xlabel("\n Tiempo transcurrido (Horas)",**hfont,fontsize=20)
# plt.ylabel("\nRadio Orbital \n Distancia Origen - Satélite (m)",**hfont,fontsize=17)
# plt.tight_layout()
# plt.xlim(0,max(t)/3600)
# plt.ylim(6200000,7500000)
# plt.legend()
# plt.show()
# plt.savefig("Grafico 4 P1.png",dpi=500)



