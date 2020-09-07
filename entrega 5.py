import xml
import xml.etree.ElementTree as ET
from numpy import zeros
import datetime as dt
from sys import argv
from time import perf_counter
import scipy as sp
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import warnings
import sys
warnings.simplefilter("ignore")



"""
===============================================================
            DEFINO LOS PARAMETROS INICIALES
===============================================================
"""



km = 1000
hora= 3600
RT = (6378.1363)*km 

            # RADIO DE LA TIERRA                – [ m ]
r = RT + (700*km)        # DISTANCIA TIERRA A SATELITE       – [ m ]
G = 6.674e-11           # CONSTANTE GRAVITACIONAL           – [ N m**2 / kg**2]
M = 5.972e24            # MASA DE LA TIERRA                 – [ kg ]
w = 7.2920150e-5         # VELOCIDAD ANGULAR DE LA TIERRA    – [ rad/s ]

atmosfera = 80*km       # ALTURA DE ATMOSFERA EN            – [ m ]

v0_x = 0         # VELOCIDAD INICIAL EN X            – [ m/s ] 
v0_y = 6820.5           # VELOCIDAD INICIAL EN Y            – [ m/s ] 
v0_z = 0                # VELOCIDAD INICIAL EN Z            – [ m/s ] 



"""
===============================================================
                  ESTABLEZCO LAS FUNCIONES
===============================================================
"""

def satelite(z,t):
    # 1. DEFINO LOS PARAMETROS Y MATRICES DE ROTACION A UTILIZAR
    θ = (w*t)
    ρ = sp.sqrt(sp.dot((z[0:3]),(z[0:3])))
    c = ((G*M)/(ρ**3))
    # c = ((G*M)/(r**3)) *******  ESTO FUE LO QUE CORREGI DE LA ENTREGA PASADA *********

    R   = sp.array([[    np.cos(θ)      ,   -np.sin(θ)    ,0       ],
                    [    np.sin(θ)      ,    np.cos(θ)    ,0       ],
                    [        0          ,        0        ,1       ]])
    
    Rp2 =2* sp.array([[   -np.sin(θ)      ,   -np.cos(θ)    ,0     ],
                      [    np.cos(θ)      ,   -np.sin(θ)    ,0     ],
                      [        0          ,        0        ,0     ]])
    
    Rpp =    sp.array([[   -np.cos(θ)      ,   np.sin(θ)    ,0     ],
                      [    -np.sin(θ)      ,  -np.cos(θ)    ,0     ],
                      [        0          ,        0        ,0     ]])

    # 2. DEFINO LOS ELEMENTOS DEL VECTOR DE ESTADO Z
    z1 = sp.array([z[0],z[1],z[2]])
    z2 = sp.array([z[3],z[4],z[5]])
    
    # 3. AQUI CREO LA TRANSFORMACION PARA EL VECTOR DE ESTADO ZP = [ZP1 , ZP2]
    zp = sp.zeros(6)
    # --- ZP1
    zp[0] = z[3] 
    zp[1] = z[4] 
    zp[2] = z[5]   
    # --- ZP2
    z2p = -(c*z1) - ( (R.T)@((w**2)*Rpp)@(z1) ) - ( (R.T)@(w*Rp2)@(z2) )
    zp[3] = z2p[0] 
    zp[4] = z2p[1] 
    zp[5] = z2p[2] 
    return zp





def satelite_mejorado(z,t):
    
    θ = (w*t)
    R   = sp.array([[    np.cos(θ)      ,   -np.sin(θ)    ,0       ],
                    [    np.sin(θ)      ,    np.cos(θ)    ,0       ],
                    [        0          ,        0        ,1       ]])
    
    Rp2 =2* sp.array([[   -np.sin(θ)      ,   -np.cos(θ)    ,0     ],
                      [    np.cos(θ)      ,   -np.sin(θ)    ,0     ],
                      [        0          ,        0        ,0     ]])
    
    Rpp =    sp.array([[   -np.cos(θ)      ,   np.sin(θ)    ,0     ],
                      [    -np.sin(θ)      ,  -np.cos(θ)    ,0     ],
                      [        0          ,        0        ,0     ]])

    x = z[0:3]
    xp = z[3:6]
    
    ρ = sp.sqrt(sp.dot(x,x))
    
    xstill = R@x
    ρnorm = (xstill)/(ρ)
    
    FG = ((G*M)/(ρ**2))*ρnorm
    FG = -((398600.440*(km**3))/(ρ**2))*ρnorm
    z2 = (xstill[2])**2
    rflat = ((xstill[0])**2) + ((xstill[1])**2)
    
    J2 = (1.7555e10)*(km**5)
    FJ2= (J2*xstill)/(ρ**7)
    FJ2[0] = (FJ2[0])*((6*z2) - (1.5*rflat))
    FJ2[1] = (FJ2[1])*((6*z2) - (1.5*rflat))
    FJ2[2] = (FJ2[2])*((3*z2) - (4.5*rflat))

    J3 = -(2.61913e11)*(km**6)
    FJ3 = sp.zeros(3)
    FJ3[0] = ((J3*xstill[0]*xstill[2]) / (ρ**9)) * (10*z2 - 7.5*rflat)
    FJ3[1] = ((J3*xstill[1]*xstill[2]) / (ρ**9)) * (10*z2 - 7.5*rflat)
    FJ3[2] = ((J3)                     / (ρ**9)) * ( (4*z2*(z2 - 3*rflat)) + (1.5*(rflat**2)))

    zp = sp.zeros(6)
    zp[0:3] = xp
    zp[3:6] = R.T@(FG + FJ2 + FJ3 - ( ((w*Rp2)@(xp)) + (((w**2)*Rpp)@(x))) )
    return zp   



def eulerint(satelite,z0,t,Nsubdivisiones):
    Nt = len(t)
    Ndim = len(sp.array(z0))
    z = zeros((Nt,Ndim))
    z[0,:] = z0
    for i in range(1,Nt):
        t_anterior = t[i-1]
        dt = (t[i] - t[i-1])/float(Nsubdivisiones)
        z_temp = z[i-1,:].copy()
        for k in range(Nsubdivisiones):
            z_temp +=  dt * satelite(z_temp, t_anterior+k*dt)
        z[i,:] = z_temp  
    return z



def utc2time(utc, ut1, EOF_datetime_format = "%Y-%m-%dT%H:%M:%S.%f"):
	t1 = dt.datetime.strptime(ut1,EOF_datetime_format)
	t2 = dt.datetime.strptime(utc,EOF_datetime_format)
	return (t2 - t1).total_seconds()


def leer_eof(fname):
	tree = ET.parse(fname)
	root = tree.getroot()

	Data_Block = root.find("Data_Block")		
	List_of_OSVs = Data_Block.find("List_of_OSVs")

	count = int(List_of_OSVs.attrib["count"])

	t = zeros(count)
	x = zeros(count)
	y = zeros(count)
	z = zeros(count)
	vx = zeros(count)
	vy = zeros(count)
	vz = zeros(count)

	set_ut1 = False
    
	for i, osv in enumerate(List_of_OSVs):
		UTC = osv.find("UTC").text[4:]
		x[i] = osv.find("X").text   #conversion de string a double es implicita
		y[i] = osv.find("Y").text
		z[i] = osv.find("Z").text
		vx[i] = osv.find("VX").text
		vy[i] = osv.find("VY").text
		vz[i] = osv.find("VZ").text
		if not set_ut1:
			ut1 = UTC
			set_ut1 = True
		t[i] = utc2time(UTC, ut1)
	return t, x, y, z, vx, vy, vz



"""
===============================================================
              ELABORO LAS VARIABLES Y GRAFICOS
===============================================================










NOTA : He comentado y descomentado secciones a lo largo de esta entrega
       asi que la gracia es que se des-comente lo que se quiera revisar
       y asi no corre todo el codigo de una! Si no se demoraria mucho. 

         AL EOF le puse max.EOF por simplicidad...






"""




# ---------- TIEMPO ------------

eofsol = leer_eof("max.EOF")
r_t = eofsol[0]
intervalo_en_segundos = r_t[-1]
t = sp.linspace(0, intervalo_en_segundos, 9361)


# ---------- ESPACIO ------------

r_x = eofsol[1]
r_y = eofsol[2]
r_z = eofsol[3]
r_vx = eofsol[4]
r_vy = eofsol[5]
r_vz = eofsol[6]

xi = r_x[0]
yi = r_y[0]
zi = r_z[0]
vxi = r_vx[0]
vyi = r_vy[0]
vzi = r_vz[0]

zi = sp.array([xi,yi,zi,vxi,vyi,vzi])



# ---------- SOLUCIONES ------------


# SOLUCION ODEINT [ o ] 
# odesol = odeint(satelite,zi,t)
# o_x = odesol[:,0]
# o_y = odesol[:,1]
# o_z = odesol[:,2]

# SOLUCION REAL [ r ]
# eofsol = leer_eof("max.EOF")
r_x = eofsol[1]
r_y = eofsol[2]
r_z = eofsol[3]

# # SOLUCION EULERINT [ e ] 
# eulersol = eulerint(satelite,zi,t,1)
# e_x = eulersol[:,0]
# e_y = eulersol[:,1]
# e_z = eulersol[:,2]


# ---------- GRAFICOS ------------



# # GRAFICO PARTE 1
# plt.figure(figsize=(12,8))
# plt.subplot(3,1,1)
# plt.title("MCOCP1-E5 MAXIMILIANO POBLETE\n TIEMPO V/S POSICION")
# plt.plot(t/hora,o_x/km    ,"salmon",linewidth=1,label="odeint")
# plt.plot(r_t/hora,r_x/km  ,"c",linewidth=1,label="real")
# plt.ylabel("X(KM)")
# plt.grid(True)

# plt.subplot(3,1,2)
# plt.plot(t/hora,o_y/km    ,"salmon",linewidth=1,label="odeint")
# plt.plot(r_t/hora,r_y/km  ,"c",linewidth=1,label="real")
# plt.ylabel("Y(KM)")
# plt.grid(True)

# plt.subplot(3,1,3)
# plt.plot(t/hora,o_z/km    ,"salmon",linewidth=1,label="odeint")
# plt.plot(r_t/hora,r_z/km  ,"c",linewidth=1,label="real")
# plt.ylabel("Z(KM)")
# plt.xlabel("Tiempo (HRS)")
# plt.grid(True)
# plt.legend()
# plt.tight_layout()
# plt.savefig("Grafico 0.png",dpi=500)
# plt.show()







# # GRAFICO PARTE 2
# odesol = odeint(satelite_mejorado,zi,t)
# o_x = odesol[:,0]
# o_y = odesol[:,1]
# o_z = odesol[:,2]
# plt.figure(figsize=(5,5))
# plt.title("MCOCP1-E5 MAXIMILIANO POBLETE\n TIEMPO V/S DERIVA ODEINT MEJORADO - SOLUCION REAL")
# delta = sp.sqrt( (r_x - o_x)**2 + (r_y - o_y)**2 + (r_z - o_z)**2 )
# plt.plot(t/hora,delta/km,"navy",label="real-odeint")
# plt.ylabel("Deriva (KM)")
# plt.xlabel("Tiempo (HRS)")
# plt.grid(True)
# plt.legend()
# plt.tight_layout()
# plt.savefig("Grafico 1.png",dpi=500)
# plt.show()
    



# # GRAFICO PARTE 3
# t1=perf_counter()
# odesol = odeint(satelite,zi,t)
# t2 = perf_counter()
# dt = t2 - t1
# print (f"\n\nEl tiempo en resolver ODEINT fue de {dt}\n\n")
# o_x = odesol[:,0]
# o_y = odesol[:,1]
# o_z = odesol[:,2]

# t1=perf_counter()
# eulersol = eulerint(satelite,zi,t,1)
# t2 = perf_counter()
# dt = t2 - t1
# print (f"\n\nEl tiempo en resolver EULERINT fue de {dt}\n\n")
# e_x = eulersol[:,0]
# e_y = eulersol[:,1]
# e_z = eulersol[:,2]

# plt.figure(figsize=(12,8))
# plt.subplot(3,1,1)
# plt.title("MCOCP1-E5 MAXIMILIANO POBLETE\n EULERINT (Nsub=1) V/S ODEINT")

# plt.plot(t/hora,o_x/km    ,"salmon",linewidth=1,label="odeint")
# plt.plot(t/hora,e_x/km  ,"yellowgreen",linewidth=1,label="eulerint")
# plt.ylabel("X(KM)")
# plt.grid(True)

# plt.subplot(3,1,2)
# plt.plot(t/hora,o_y/km    ,"salmon",linewidth=1,label="odeint")
# plt.plot(t/hora,e_y/km  ,"yellowgreen",linewidth=1,label="eulerint")
# plt.ylabel("Y(KM)")
# plt.grid(True)

# plt.subplot(3,1,3)
# plt.plot(t/hora,o_z/km    ,"salmon",linewidth=1,label="odeint")
# plt.plot(t/hora,e_z/km  ,"yellowgreen",linewidth=1,label="eulerint (Nsub=1)")
# plt.ylabel("Z(KM)")
# plt.xlabel("Tiempo (HRS)")
# plt.grid(True)
# plt.legend()
# plt.tight_layout()
# plt.savefig("Grafico 2.png",dpi=500)
# plt.show()






# # GRAFICO PARTE 4


# plt.figure(figsize=(5,5))
# plt.title("MCOCP1-E5 MAXIMILIANO POBLETE\n TIEMPO V/S DERIVA EULERINT (Nsub=1) - ODEINT")
# delta = sp.sqrt( (e_x - o_x)**2 + (e_y - o_y)**2 + (e_z - o_z)**2 )
# plt.plot(t/hora,delta/km,"navy",label="eulerint (Nsub=1) - odeint ")
# plt.ylabel("Deriva (KM)")
# plt.xlabel("Tiempo (HRS)")
# plt.grid(True)
# plt.legend()
# plt.tight_layout()
# plt.savefig("Grafico 3.png",dpi=500)
# plt.show()
# print (f"La deriva al final, entre ODEINT y EULERINT (Nsub=1) , es de {delta[-1]/km} km")






# CALCULO DE ERROR PARTE 1

# Nsubs = 2

# eulersol = eulerint(satelite,zi,t,Nsubs)
# e_xf = eulersol[:,0][-1]
# e_yf = eulersol[:,1][-1]
# e_zf = eulersol[:,2][-1]

# r_xf = eofsol[1][-1]
# r_yf = eofsol[2][-1]
# r_zf = eofsol[3][-1]

# e_delta = sp.sqrt( (e_xf)**2 + (e_yf)**2 + (e_zf)**2 )
# r_delta = sp.sqrt( (r_xf)**2 + (r_yf)**2 + (r_zf)**2 )

# # error = ((abs(e_delta - r_delta))/(r_delta))*100
# error = 100

# while error > 1:
#     t1=perf_counter()
#     eulersol = eulerint(satelite,zi,t,Nsubs)
#     t2 = perf_counter()
#     dt = t2 - t1
#     e_xf = eulersol[:,0][-1]
#     e_yf = eulersol[:,1][-1]
#     e_zf = eulersol[:,2][-1]
#     e_delta = sp.sqrt( (e_xf)**2 + (e_yf)**2 + (e_zf)**2 )    
#     error = ((abs(e_delta - r_delta))/(r_delta))*100
#     print(f"Para {Nsubs} el tiempo es {dt} segundos y el error es de {error}%")
#     Nsubs += 1
# print (f"Finalmente, Para {Nsubs} el tiempo es {dt} segundos y el error es de {error}%")








# CALCULO DE ERROR PARTE 2

# Nsubs = 10
# eulersol = eulerint(satelite,zi,t,Nsubs)
# e_xf = eulersol[:,0][-1]
# e_yf = eulersol[:,1][-1]
# e_zf = eulersol[:,2][-1]
# e_x = eulersol[:,0]
# e_y = eulersol[:,1]
# e_z = eulersol[:,2]

# r_xf = eofsol[1][-1]
# r_yf = eofsol[2][-1]
# r_zf = eofsol[3][-1]

# distancia = sp.sqrt( (e_xf-r_xf)**2 + (e_yf-r_yf)**2 +  (e_zf-r_zf)**2 )
# e_delta = sp.sqrt( (e_xf)**2 + (e_yf)**2 + (e_zf)**2 )
# r_delta = sp.sqrt( (r_xf)**2 + (r_yf)**2 + (r_zf)**2 )
# error = (distancia / r_delta)*100
# print (error)
# while error > 1:
#     eulersol = eulerint(satelite,zi,t,Nsubs)
#     e_xf = eulersol[:,0][-1]
#     e_yf = eulersol[:,1][-1]
#     e_zf = eulersol[:,2][-1]
#     distancia = sp.sqrt( (e_xf-r_xf)**2 + (e_yf-r_yf)**2 +  (e_zf-r_zf)**2 )
#     e_delta = sp.sqrt( (e_xf)**2 + (e_yf)**2 + (e_zf)**2 )
#     r_delta = sp.sqrt( (r_xf)**2 + (r_yf)**2 + (r_zf)**2 )
#     error = (distancia / r_delta)*100
#     print(f"Para {Nsubs} el tiempo es {dt} segundos y el error es de {error}%")
#     Nsubs += 1









# # # GRAFICO PARTE 5
# Nsubs = 200
# t1=perf_counter()
# eulersol = eulerint(satelite,zi,t,Nsubs)
# t2 = perf_counter()
# dt = t2 - t1
# print (f"\n\nEl tiempo en resolver EULERINT con Nsubs = {Nsubs} fue de {dt}\n\n")
# e_x = eulersol[:,0]
# e_y = eulersol[:,1]
# e_z = eulersol[:,2]

# plt.figure(figsize=(5,5))
# plt.title(f"MCOCP1-E5 MAXIMILIANO POBLETE\n TIEMPO V/S DERIVA EULERINT (Nsub={Nsubs}) - SOLUCION REAL")
# delta = sp.sqrt( (e_x - r_x)**2 + (e_y - r_y)**2 + (e_z - r_z)**2 )
# plt.plot(t/hora,delta/km,"navy",label="eulerint (Nsub=200) - real ")
# plt.ylabel("Deriva (KM)")
# plt.xlabel("Tiempo (HRS)")
# plt.grid(True)
# plt.legend()
# plt.tight_layout()
# plt.savefig("Grafico 4.png",dpi=500)
# plt.show()
# print (f"La deriva al final, entre REAL y EULERINT (Nsub={Nsubs}) , es de {delta[-1]/km} km")








# # # GRAFICO PARTE 6

# plt.figure(figsize=(12,8))
# plt.subplot(3,1,1)
# plt.title("MCOCP1-E5 MAXIMILIANO POBLETE\n EULERINT (Nsub=200) V/S REAL")

# plt.plot(t/hora,r_x/km    ,"salmon",linewidth=1,label="real")
# plt.plot(t/hora,e_x/km  ,"yellowgreen",linewidth=1,label="eulerint")
# plt.ylabel("X(KM)")
# plt.grid(True)

# plt.subplot(3,1,2)
# plt.plot(t/hora,r_y/km    ,"salmon",linewidth=1,label="real")
# plt.plot(t/hora,e_y/km  ,"yellowgreen",linewidth=1,label="eulerint")
# plt.ylabel("Y(KM)")
# plt.grid(True)

# plt.subplot(3,1,3)
# plt.plot(t/hora,r_z/km    ,"salmon",linewidth=1,label="real")
# plt.plot(t/hora,e_z/km  ,"yellowgreen",linewidth=1,label="eulerint (Nsub=200)")
# plt.ylabel("Z(KM)")
# plt.xlabel("Tiempo (HRS)")
# plt.grid(True)
# plt.legend()
# plt.tight_layout()
# plt.savefig("Grafico 5.png",dpi=500)
# plt.show()











# # ---------- SOLUCIONES PARA PARTE 4 ------------

# t0total = perf_counter()


# # ---------- TIEMPO ------------
# eofsol = leer_eof("max.EOF")
# r_t = eofsol[0]
# intervalo_en_segundos = r_t[-1]
# t = sp.linspace(0, intervalo_en_segundos, 9361)


# # ---------- ESPACIO ------------
# r_x = eofsol[1]
# r_y = eofsol[2]
# r_z = eofsol[3]
# r_vx = eofsol[4]
# r_vy = eofsol[5]
# r_vz = eofsol[6]

# xi = r_x[0]
# yi = r_y[0]
# zi = r_z[0]
# vxi = r_vx[0]
# vyi = r_vy[0]
# vzi = r_vz[0]

# zi = sp.array([xi,yi,zi,vxi,vyi,vzi])



# #SOLUCION ODEINT MEJORADA [ o ] 
# t1=perf_counter()
# odesol = odeint(satelite_mejorado,zi,t)
# t2=perf_counter()
# o_x = odesol[:,0]
# o_y = odesol[:,1]
# o_z = odesol[:,2]

# # SOLUCION REAL [ r ]
# eofsol = leer_eof("max.EOF")
# r_x = eofsol[1]
# r_y = eofsol[2]
# r_z = eofsol[3]

# # ---------- GRAFICOS ------------

# # GRAFICO PARTE 7

# plt.figure(figsize=(12,8))
# plt.subplot(3,1,1)
# plt.title("MCOCP1-E5 MAXIMILIANO POBLETE\n TIEMPO V/S POSICION")
# plt.plot(t/hora,o_x/km    ,"salmon",linewidth=1,label="odeint mejorada")
# plt.plot(r_t/hora,r_x/km  ,"c",linewidth=1,label="real")
# plt.ylabel("X(KM)")
# plt.grid(True)

# plt.subplot(3,1,2)
# plt.plot(t/hora,o_y/km    ,"salmon",linewidth=1,label="odeint mejorada")
# plt.plot(r_t/hora,r_y/km  ,"c",linewidth=1,label="real")
# plt.ylabel("Y(KM)")
# plt.grid(True)

# plt.subplot(3,1,3)
# plt.plot(t/hora,o_z/km    ,"salmon",linewidth=1,label="odeint mejorada")
# plt.plot(r_t/hora,r_z/km  ,"c",linewidth=1,label="real")
# plt.ylabel("Z(KM)")
# plt.xlabel("Tiempo (HRS)")
# plt.grid(True)
# plt.legend()
# plt.tight_layout()
# plt.savefig("Grafico 6.png",dpi=500)
# plt.show()


# # GRAFICO PARTE 8
# plt.figure(figsize=(5,5))
# plt.title("MCOCP1-E5 MAXIMILIANO POBLETE\n TIEMPO V/S DERIVA ODEINT MEJORADO - SOLUCION REAL")
# delta = sp.sqrt( (r_x - o_x)**2 + (r_y - o_y)**2 + (r_z - o_z)**2 )
# plt.plot(t/hora,delta/km,"navy",label="real-odeint")
# plt.ylabel("Deriva (KM)")
# plt.xlabel("Tiempo (HRS)")
# plt.grid(True)
# plt.legend()
# plt.tight_layout()
# plt.savefig("Grafico 7.png",dpi=500)
# plt.show()
    
# tftotal = perf_counter()
# print("El tiempo transucrrido total fue de",(tftotal-t0total))
# print("El tiempo transucrrido para la solucion de odeint mejorado fue de",(t2-t1))