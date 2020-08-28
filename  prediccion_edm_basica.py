import scipy as sp
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import warnings
import sys
warnings.simplefilter("ignore")


"""


███████╗██╗░░░░░      ░█████╗░░█████╗░██████╗░██╗░██████╗░░█████╗░    
██╔════╝██║░░░░░      ██╔══██╗██╔══██╗██╔══██╗██║██╔════╝░██╔══██╗    
█████╗░░██║░░░░░      ██║░░╚═╝██║░░██║██║░░██║██║██║░░██╗░██║░░██║    
██╔══╝░░██║░░░░░      ██║░░██╗██║░░██║██║░░██║██║██║░░╚██╗██║░░██║    
███████╗███████╗      ╚█████╔╝╚█████╔╝██████╔╝██║╚██████╔╝╚█████╔╝    
╚══════╝╚══════╝      ░╚════╝░░╚════╝░╚═════╝░╚═╝░╚═════╝░░╚════╝░    

██████╗░███████╗      ███████╗░██████╗████████╗░█████╗░    
██╔══██╗██╔════╝      ██╔════╝██╔════╝╚══██╔══╝██╔══██╗    
██║░░██║█████╗░░      █████╗░░╚█████╗░░░░██║░░░███████║    
██║░░██║██╔══╝░░      ██╔══╝░░░╚═══██╗░░░██║░░░██╔══██║    
██████╔╝███████╗      ███████╗██████╔╝░░░██║░░░██║░░██║    
╚═════╝░╚══════╝      ╚══════╝╚═════╝░░░░╚═╝░░░╚═╝░░╚═╝    

███████╗███╗░░██╗████████╗██████╗░███████╗░██████╗░░█████╗░      ██████╗░  
██╔════╝████╗░██║╚══██╔══╝██╔══██╗██╔════╝██╔════╝░██╔══██╗      ╚════██╗  
█████╗░░██╔██╗██║░░░██║░░░██████╔╝█████╗░░██║░░██╗░███████║      ░█████╔╝  
██╔══╝░░██║╚████║░░░██║░░░██╔══██╗██╔══╝░░██║░░╚██╗██╔══██║      ░╚═══██╗  
███████╗██║░╚███║░░░██║░░░██║░░██║███████╗╚██████╔╝██║░░██║      ██████╔╝  
╚══════╝╚═╝░░╚══╝░░░╚═╝░░░╚═╝░░╚═╝╚══════╝░╚═════╝░╚═╝░░╚═╝      ╚═════╝░  

███████╗███╗░░░███╗██████╗░██╗███████╗███████╗░█████╗░      ███████╗███╗░░██╗    
██╔════╝████╗░████║██╔══██╗██║██╔════╝╚════██║██╔══██╗      ██╔════╝████╗░██║    
█████╗░░██╔████╔██║██████╔╝██║█████╗░░░░███╔═╝███████║      █████╗░░██╔██╗██║    
██╔══╝░░██║╚██╔╝██║██╔═══╝░██║██╔══╝░░██╔══╝░░██╔══██║      ██╔══╝░░██║╚████║    
███████╗██║░╚═╝░██║██║░░░░░██║███████╗███████╗██║░░██║      ███████╗██║░╚███║    
╚══════╝╚═╝░░░░░╚═╝╚═╝░░░░░╚═╝╚══════╝╚══════╝╚═╝░░╚═╝      ╚══════╝╚═╝░░╚══╝    

██╗░░░░░░█████╗░      ██╗░░░░░██╗███╗░░██╗███████╗░█████╗░    
██║░░░░░██╔══██╗      ██║░░░░░██║████╗░██║██╔════╝██╔══██╗    
██║░░░░░███████║      ██║░░░░░██║██╔██╗██║█████╗░░███████║    
██║░░░░░██╔══██║      ██║░░░░░██║██║╚████║██╔══╝░░██╔══██║    
███████╗██║░░██║      ███████╗██║██║░╚███║███████╗██║░░██║    
╚══════╝╚═╝░░╚═╝      ╚══════╝╚═╝╚═╝░░╚══╝╚══════╝╚═╝░░╚═╝    

██████╗░███████╗░█████╗░
╚════██╗╚════██║██╔══██╗
░░███╔═╝░░░░██╔╝██║░░██║
██╔══╝░░░░░██╔╝░██║░░██║
███████╗░░██╔╝░░╚█████╔╝
╚══════╝░░╚═╝░░░░╚════╝░

"""












"""

███████╗███╗░░██╗████████╗██████╗░███████╗░██████╗░░█████╗░      ██████╗░
██╔════╝████╗░██║╚══██╔══╝██╔══██╗██╔════╝██╔════╝░██╔══██╗      ╚════██╗
█████╗░░██╔██╗██║░░░██║░░░██████╔╝█████╗░░██║░░██╗░███████║      ░░███╔═╝
██╔══╝░░██║╚████║░░░██║░░░██╔══██╗██╔══╝░░██║░░╚██╗██╔══██║      ██╔══╝░░
███████╗██║░╚███║░░░██║░░░██║░░██║███████╗╚██████╔╝██║░░██║      ███████╗
╚══════╝╚═╝░░╚══╝░░░╚═╝░░░╚═╝░░╚═╝╚══════╝░╚═════╝░╚═╝░░╚═╝      ╚══════╝

"""


"""
===============================================================
            DEFINO LOS PARAMETROS INICIALES
===============================================================
"""
km = 1000
hora=3600
RT = 6371*km             # RADIO DE LA TIERRA                – [ m ]
r = RT + (700*km)        # DISTANCIA TIERRA A SATELITE       – [ m ]
G = 6.674e-11           # CONSTANTE GRAVITACIONAL           – [ N m**2 / kg**2]
M = 5.972e24            # MASA DE LA TIERRA                 – [ kg ]
w = 7.272205e-5         # VELOCIDAD ANGULAR DE LA TIERRA    – [ rad/s ]

atmosfera = 80*km       # ALTURA DE ATMOSFERA EN            – [ m ]

v0_x = 0         # VELOCIDAD INICIAL EN X            – [ m/s ] 
v0_y = 6820.5           # VELOCIDAD INICIAL EN Y            – [ m/s ] 
v0_z = 0                # VELOCIDAD INICIAL EN Z            – [ m/s ] 



"""
===============================================================
         ESTABLEZCO LA FUNCION PARA LAS ORBITAS
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




"""
===============================================================
            IMPORTO LAS LIBRERIAS PARA GRAFICAR
===============================================================
"""


#****************************************************************************************************
#****************************************************************************************************
#****************************************************************************************************


plt.figure(1,figsize=(8,9))

"----  GRAFICO LA ORBITA  ----"

t = sp.linspace(0,3.2*hora,1001)
z0  = sp.array([r,0,0,v0_x,v0_y,v0_z])
sol = odeint(satelite,z0,t)
x = sol[:,0]
y = sol[:,1]
z = sol[:,2]

plt.plot(x,y,"orange",linewidth=1.5,label=f"V = {v0_y} m/s ")
plt.plot(r,0,"bo",markersize=5,label="Inicio")



"----  GRAFICO LA TIERRA  ----"
α = sp.linspace(0,360,360)
plt.plot(RT*np.cos(np.radians(α)),RT*np.sin(np.radians(α)),"c")
Rplot = RT
for i in range(50):
    Rplot=Rplot-5700*i
    plt.plot(Rplot*np.cos(np.radians(α)),Rplot*np.sin(np.radians(α)),"c",alpha=0.2)




"----  GRAFICO LA ATMOSFERA  ----"
plt.plot((RT+atmosfera)*np.cos(np.radians(α)),(RT+atmosfera)*np.sin(np.radians(α)),"--r",linewidth=1,label="Atmosfera")




"----   EJES DE LA TIERRA     ----"

plt.vlines(0,-RT,RT)
plt.hlines(0,-RT,RT)




"----  PARAMETROS EXTRAS PARA GRAFICAR   ----"



hfont = {'fontname':'Avenir'}
plt.title("MCOCP1-E2 MAXIMILIANO POBLETE\n ORBITA SATELITAL",**hfont,fontsize=20)


plt.xticks(**hfont,rotation=90,fontsize=15)
plt.xticks(**hfont,fontsize=16)
plt.yticks(**hfont,fontsize=16)


plt.xlabel("Posicion en X (m)",**hfont,fontsize=17)
plt.ylabel("Posicion en Y (m)",**hfont,fontsize=17)

plt.ticklabel_format(useOffset=False, style='plain')

plt.grid(True)

plt.tight_layout()

plt.legend()


# plt.savefig("Grafico 0 P1.png",dpi=500)

# plt.show()
plt.close()



#****************************************************************************************************
#****************************************************************************************************
#****************************************************************************************************



plt.figure(2,figsize=(8,9))



radioactual=[]
for i in range(len(x)):
    radioactual.append(np.sqrt (x[i]**2 + y[i]**2 + z[i]**2))
plt.plot(t,radioactual)
plt.hlines(  RT ,0,max(t),"g",linewidth=5)
plt.hlines(  RT+atmosfera ,0,max(t),"r","--")






plt.title("MCOCP1-E2 MAXIMILIANO POBLETE\n TIEMPO V/S RADIO ORBITAL",**hfont,fontsize=20)
plt.ticklabel_format(useOffset=False, style='plain')
plt.xticks(**hfont,rotation=90,fontsize=15)
plt.xticks(**hfont,fontsize=16)
plt.yticks(**hfont,fontsize=16)
plt.grid(True)

plt.xlabel("Tiempo transcurrido (Horas)",**hfont,fontsize=17)
plt.ylabel("Radio Orbital – Distancia Origen - Satélite (m)",**hfont,fontsize=17)
# plt.show()
plt.close()


#****************************************************************************************************
#****************************************************************************************************
#****************************************************************************************************

"""

███████╗███╗░░██╗████████╗██████╗░███████╗░██████╗░░█████╗░      ██████╗░
██╔════╝████╗░██║╚══██╔══╝██╔══██╗██╔════╝██╔════╝░██╔══██╗      ╚════██╗
█████╗░░██╔██╗██║░░░██║░░░██████╔╝█████╗░░██║░░██╗░███████║      ░█████╔╝
██╔══╝░░██║╚████║░░░██║░░░██╔══██╗██╔══╝░░██║░░╚██╗██╔══██║      ░╚═══██╗
███████╗██║░╚███║░░░██║░░░██║░░██║███████╗╚██████╔╝██║░░██║      ██████╔╝
╚══════╝╚═╝░░╚══╝░░░╚═╝░░░╚═╝░░╚═╝╚══════╝░╚═════╝░╚═╝░░╚═╝      ╚═════╝░


"""





import datetime as dt





"""

░█████╗░░█████╗░██████╗░██╗░██████╗░░█████╗░      ███╗░░░███╗██╗░█████╗░
██╔══██╗██╔══██╗██╔══██╗██║██╔════╝░██╔══██╗      ████╗░████║██║██╔══██╗
██║░░╚═╝██║░░██║██║░░██║██║██║░░██╗░██║░░██║      ██╔████╔██║██║██║░░██║
██║░░██╗██║░░██║██║░░██║██║██║░░╚██╗██║░░██║      ██║╚██╔╝██║██║██║░░██║
╚█████╔╝╚█████╔╝██████╔╝██║╚██████╔╝╚█████╔╝      ██║░╚═╝░██║██║╚█████╔╝
░╚════╝░░╚════╝░╚═════╝░╚═╝░╚═════╝░░╚════╝░      ╚═╝░░░░░╚═╝╚═╝░╚════╝░

"""


"""

----------------------------------------------------------------------
  ----------   DATOS EOF PARA MAXIMILIANO POBLETE ----------------
----------------------------------------------------------------------



PRIMER COMPONENTE

  <List_of_OSVs count="9361">
    <OSV>
      <TAI>TAI=2020-07-25T23:00:19.000000</TAI>
      <UTC>UTC=2020-07-25T22:59:42.000000</UTC>
      <UT1>UT1=2020-07-25T22:59:41.787032</UT1>
      <Absolute_Orbit>+22634</Absolute_Orbit>
      <X unit="m">-394791.420488</X>
      <Y unit="m">-2222608.647191</Y>
      <Z unit="m">6695831.133945</Z>
      <VX unit="m/s">-2343.275744</VX>
      <VY unit="m/s">6887.111751</VY>
      <VZ unit="m/s">2143.280850</VZ>
      <Quality>NOMINAL</Quality>
    </OSV>



ULTIMO COMPONENTE

<OSV>
      <TAI>TAI=2020-07-27T01:00:19.000000</TAI>
      <UTC>UTC=2020-07-27T00:59:42.000000</UTC>
      <UT1>UT1=2020-07-27T00:59:41.787584</UT1>
      <Absolute_Orbit>+22650</Absolute_Orbit>
      <X unit="m">-1818680.097930</X>
      <Y unit="m">-6838669.773573</Y>
      <Z unit="m">72865.206454</Z>
      <VX unit="m/s">-1508.337828</VX>
      <VY unit="m/s">489.698215</VY>
      <VZ unit="m/s">7430.040379</VZ>
      <Quality>NOMINAL</Quality>
    </OSV>

"""

utc_EOF_format = "%Y-%m-%dT%H:%M:%S.%f"
t1 = dt.datetime.strptime("2020-07-25T22:59:42.000000",utc_EOF_format)
t2 = dt.datetime.strptime("2020-07-27T00:59:42.000000",utc_EOF_format)

intervalo = t2-t1
intervalo_en_segundos = intervalo.total_seconds()


xi = -394791.420488
yi = -2222608.647191
zi = 6695831.133945

vxi = -2343.275744
vyi = 6887.111751
vzi = 2143.280850


xf = -1818680.097930
yf = -6838669.773573
zf = 72865.206454

vxf = -1508.337828
vyf = 489.698215
vzf = 7430.040379


q = sp.linspace(0,  intervalo_en_segundos  ,93600)


zi = sp.array([ xi,
                yi,
                zi,
                vxi,
                vyi,
                vzi ])

zf = sp.array([ xf,
                yf,
                zf,
                vxf,
                vyf,
                vzf])


sol = odeint(satelite, zi, q)

pos_final = zf - sol[-1]

dx = pos_final[0]
dy = pos_final[1]
dz = pos_final[2]

print (np.sqrt(dx**2 + dy**2 + dz**2))







# """

# ░█████╗░░█████╗░██████╗░██╗░██████╗░░█████╗░        ██████╗░██████╗░░█████╗░███████╗███████╗
# ██╔══██╗██╔══██╗██╔══██╗██║██╔════╝░██╔══██╗        ██╔══██╗██╔══██╗██╔══██╗██╔════╝██╔════╝
# ██║░░╚═╝██║░░██║██║░░██║██║██║░░██╗░██║░░██║        ██████╔╝██████╔╝██║░░██║█████╗░░█████╗░░
# ██║░░██╗██║░░██║██║░░██║██║██║░░╚██╗██║░░██║        ██╔═══╝░██╔══██╗██║░░██║██╔══╝░░██╔══╝░░
# ╚█████╔╝╚█████╔╝██████╔╝██║╚██████╔╝╚█████╔╝        ██║░░░░░██║░░██║╚█████╔╝██║░░░░░███████╗
# ░╚════╝░░╚════╝░╚═════╝░╚═╝░╚═════╝░░╚════╝░        ╚═╝░░░░░╚═╝░░╚═╝░╚════╝░╚═╝░░░░░╚══════╝


# """


# """

# ----------------------------------------------------------------------
#   ----------   DATOS EOF DE CANVAS           ----------------
# ----------------------------------------------------------------------


# Vector de estado inicial:

#   TAI=2018-08-14T23:00:19.000000
#   UTC=2018-08-14T22:59:42.000000
#   UT1=2018-08-14T22:59:42.065788
#   +23248
#   2083293.582654
#   -6380690.028717
#   -2250417.463178
#   -805.237768
#   -2737.127661
#   7035.393528
  


# Vector de estado final:

#   TAI=2018-08-16T01:00:19.000000
#   UTC=2018-08-16T00:59:42.000000
#   UT1=2018-08-16T00:59:42.065173
#   +23264
#   1010095.475188
#   -123073.229951
#   -7008913.872221
#   -1887.279838
#   -7324.093741
#   -143.308102
  
# """











# t3 = dt.datetime.strptime("2018-08-14T22:59:42.000000",utc_EOF_format)
# t4 = dt.datetime.strptime("2018-08-16T00:59:42.000000",utc_EOF_format)

# intervalo = t4-t3
# intervalo_en_segundos_profesor = intervalo.total_seconds()


# xi=  2083293.582654
# yi=  -6380690.028717
# zi=  -2250417.463178
# vxi=  -805.237768
# vyi= -2737.127661
# vzi=  7035.393528

# xf =  1010095.475188
# yf =  -123073.229951
# zf = -7008913.872221
# vxf =  -1887.279838
# vyf =  -7324.093741
# vzf =  -143.308102

# p = sp.linspace(0,intervalo_en_segundos_profesor,93600)

# zi = sp.array([ xi,
#                 yi,
#                 zi,
#                 vxi,
#                 vyi,
#                 vzi ])

# zf = sp.array([ xf,
#                 yf,
#                 zf,
#                 vxf,
#                 vyf,
#                 vzf])


# sol_profesor = odeint(satelite, zi, p)

# pos_final_profesor = zf - sol_profesor[-1]

# dx = pos_final_profesor[0]
# dy = pos_final_profesor[1]
# dz = pos_final_profesor[2]

# print (np.sqrt(dx**2 + dy**2 + dz**2))































      