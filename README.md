# MCOC2020-P1

### E1 - Integración de ecuaciones diferenciales (Lanzamiento de Proyectil)

![alt text](https://github.com/maxipoblete/MCOC2020-P1/blob/master/Proyectiles.png )

### E2 -  Primeras predicciones con la EDM básica del satélite
Entregado en Canvas

### E3 -  I/O de vectores de estado y predicciones usando la EDM básica

1849958.4027598512
### E4 -  Estudio de convergencia Método de Euler
Entregado en Canvas

### E5 -  I/O de vectores de estado y predicciones usando la EDM básica

#### 1. Grafíque, como arriba, la posición (x,y,z) en el tiempo del vector de estado de Sentinel 1A/B que le tocó. Para esto, descargue y utilice la función leer_eof.py (Enlaces a un sitio externo.) para poder trabajar con los archivos EOF.

Por motivos personales estaré utilizando mi propia sintaxis de codigo, según las 3 soluciones que se deben analizar. La solucion real EOF ( r ), La solucion odeint ( o ), y la Solucion Eulerint ( e ); todas estas como prefijo a las variables. De esta manera, el grafico solicitado es el siguiente. La curva de color azul indica la solucion real, y la curva de color salmon indica la solucion odeint básica (sin mejorar).

![alt text](   https://github.com/maxipoblete/MCOC2020-P1/blob/master/Grafico%200.png   )



#### 2. Usando la condición inicial (primer OSV) de su archivo, compare la solución entre odeint y eulerint. Use Nsubdiviciones=1. Grafíque la deriva en el tiempo como arriba ¿Cuánto deriva eulerint de odeint en este caso al final del tiempo? (Esta pregunta solo compara algoritmos, no se usa más que la condición inicial del archivo EOF). ¿Cuanto se demora odeint y eulerint respectivamente en producir los resultados?

Para realizar esto, se utiliza el modulo time y la funcion perf_counter(). De esta forma, se obtiene que el tiempo de resolucion es el siguiente, con Nsubdivisiones=1:<br>
El tiempo en resolver ODEINT fue de 0.2889721940000527 s <br>
El tiempo en resolver EULERINT fue de 0.891989900000226 s <br>
Luego, se obtuvo que la deriva al final del tiempo, entre EULERINT y ODEINT, medido en kilometros fue de 20126.86813146657 km <br>
<br>
A continuacion se muestran los graficos obtenidos.
![alt text](   https://github.com/maxipoblete/MCOC2020-P1/blob/master/Grafico%202.png  )
![alt text](   https://github.com/maxipoblete/MCOC2020-P1/blob/master/Grafico%203.png  )

#### 3. ¿Cuantas subdivisiones hay que usar para que la predicción con eulerint al final del tiempo esté en menos de un 1% de error? ⚠️[SE ASUME QUE SE DEBE CALCULAR ENTRE LA REAL Y LA SOLUCION EULERINT] Grafique la deriva en el tiempo. Comente con respecto del tiempo de ejecución de eulerint ahora. 

Cada vez que se incrementa en un paso una subdivision, el tiempo de ejecucion incrementa considerablemente. Para un N=200 se demoró 177.32198137100022 segundos y la deriva final fue de 3662.7733398098885 km. El error fue calculado como el error entre la distancia del punto final para los valores reales y la solucion de euler int, con la distancia del radio desde el origen hasta el punto final de la solucion real. Con esto se obtiene un 50% de error hacia el final del resultado. Lo que se puede apreciar en los siguientes graficos:

![alt text](   https://github.com/maxipoblete/MCOC2020-P1/blob/master/Grafico%2044.png  )
![alt text](   https://github.com/maxipoblete/MCOC2020-P1/blob/master/Grafico%2055.png  )

Se ve que para aproximadamente la mitad del tiempo la deriva es baja pero aumenta explosivamente a partir de 10 horas. Con el calculo de error, se tendrán valores oscilantes a medida que aumentamos la cantidad de Nsubdivisiones hasta un punto en que el grafico de la prediccion se apegue al grafico de la solucion real completamente. Esto puede tardar muchisimo rato. Para valores de N=1500 , el error fue de aproximadamente 11% y se demoro como 30 minutos. (No alcance a graficarlo ya que estaba probando). Con valores de N=3000 a N=5000 subdivisiones, debería ser suficiente para obtener errores menores a lo solicitado. (Lamento que mi Pc no sea tan potente para haber probado esto bien, en el codigo se ve que al menos lo intente).


#### 4. Implemente las correcciones J2 y J3. Repita los gráficos de arriba para su caso particular. ¿Cuánta deriva incurre al agregar las correcciones J2 y J3? ¿Cuanto se demora su código en correr?

En esta parte, se implemento la con los términos J2 y J3. Es decir, se comparará entre la Solucion Real y la solucion Odeint Mejorada.<br>
<br> 
En primer lugar, se puede obtener que la demora en correr el codigo para la resolucion de Odeint Mejorado es de: 0.3466004299989436 s<br>
Si se considera el tiempo total desde que la lectura del archivo hasta a la ultima linea , se obtienen: 3.8667460739998205 s <br>
<br>
A continuacion se presenta el grafico de las posiciones y la deriva entre Odeint Mejorado y la Solucion Real:


![alt text](   https://github.com/maxipoblete/MCOC2020-P1/blob/master/Grafico%206.png  )
![alt text](   https://github.com/maxipoblete/MCOC2020-P1/blob/master/Grafico%207.png  )

Se obtuvo una mejora considerable! La deriva al final resulto ser de 4 km. En el grafico de las posiciones , se puede apreciar mucho mejor si se hace un zoom hacia la parte final:

![alt text](   https://github.com/maxipoblete/MCOC2020-P1/blob/master/Grafico%208.png  )

<br>
<br>
<br>
<br>

## ✅💥✅💥 E7 -  Entrega Final ✅💥 ✅💥

#### Reporte de Ajustes realizados

En esta entrega se han agregado 8 terminos zonales y teserales, utilizando el procedimiento descrito por el profesor en el video que subio de SymPy:
1. Se tiene una funcion ***legendre(m,n,x)*** que recibe los coeficientes m, n y devuelve las funciones de legendre en funcion de la variable x. En este caso la variable de entrada corresponde a **sinθ=z/r** donde r corresponde a la distancia desde el origen hasta el satelite (esfericas) r=sqrt(x²+y²+z²). La idea es definir en SymPy las variables que se utilizaran en la fórmula del potencial "u", el cual se deriva para obtener la fuerza Fg que será incluida en la ecuación a solucionar por la predicción *odeint*. Estas variables son:


```
x = S('x') -----> Coordenada x
y = S('y') -----> Coordenada y
z = S('z') -----> Coordenada z
μ= S('μ')  -----> 398600441500000.000000 m^3/s^2 : G*(masa de la tierra)
J2=S('J2')  -----> Terminos Zonales JX...
J3=S('J3')
J4=S('J4')
J5=S('J5')
J6=S('J6')
J7=S('J7')
J8 =S('J8')
R=S('R')    -----> 6378136.3  m  : Radio de la tierra
```

Luego de tener las variables listas en formato simbolico, basta con definir el potencial "u" según la fórmula del modelo de Wikipedia y derivar parcialmente en cada componente:

```
sin0 = z/r
φ = atan2(y,x)     


u = -μ/r*(1+
          (  J2 * legendre(0,2,sin0) / (r/R)**2 ) +   (  (legendre(1,2,sin0))*( (C[2,1])*cos(1*φ) + (S[2,1])*sin(1*φ)) / ((r/R)**2) ) +  ( (legendre(2,2,sin0))*( (C[2,2])*cos(2*φ) + (S[2,2])*sin(2*φ))/((r/R)**2) )
          + (J3 * legendre(0,3,sin0) / (r/R)**3) +    (  (legendre(1,3,sin0))*( (C[3,1])*cos(1*φ) + (S[3,1])*sin(1*φ)) / ((r/R)**3) ) + (  (legendre(2,3,sin0))*( (C[3,2])*cos(2*φ) + (S[3,2])*sin(2*φ)) / ((r/R)**3) ) + (  (legendre(3,3,sin0))*( (C[3,3])*cos(3*φ) + (S[3,3])*sin(3*φ)) / ((r/R)**3) )
          + (J4 * legendre(0,4,sin0) / (r/R)**4) +    (  (legendre(1,4,sin0))*( (C[4,1])*cos(1*φ) + (S[4,1])*sin(1*φ)) / ((r/R)**4) ) + (  (legendre(2,4,sin0))*( (C[4,2])*cos(2*φ) + (S[4,2])*sin(2*φ)) / ((r/R)**4) ) + (  (legendre(3,4,sin0))*( (C[4,3])*cos(3*φ) + (S[4,3])*sin(3*φ)) / ((r/R)**4) ) +  (  (legendre(4,4,sin0))*( (C[4,4])*cos(4*φ) + (S[4,4])*sin(4*φ)) / ((r/R)**4) )  
          + (J5 * legendre(0,5,sin0) / (r/R)**5) +    (  (legendre(1,5,sin0))*( (C[5,1])*cos(1*φ) + (S[5,1])*sin(1*φ)) / ((r/R)**5) ) + (  (legendre(2,5,sin0))*( (C[5,2])*cos(2*φ) + (S[5,2])*sin(2*φ)) / ((r/R)**5) ) + (  (legendre(3,5,sin0))*( (C[5,3])*cos(3*φ) + (S[5,3])*sin(3*φ)) / ((r/R)**5) ) +  (  (legendre(4,5,sin0))*( (C[5,4])*cos(4*φ) + (S[5,4])*sin(4*φ)) / ((r/R)**5) )+  (  (legendre(5,5,sin0))*( (C[5,5])*cos(5*φ) + (S[5,5])*sin(5*φ)) / ((r/R)**5) ) 
          + (J6 * legendre(0,6,sin0) / (r/R)**6) +    (  (legendre(1,6,sin0))*( (C[6,1])*cos(1*φ) + (S[6,1])*sin(1*φ)) / ((r/R)**6) ) + (  (legendre(2,6,sin0))*( (C[6,2])*cos(2*φ) + (S[6,2])*sin(2*φ)) / ((r/R)**6) ) + (  (legendre(3,6,sin0))*( (C[6,3])*cos(3*φ) + (S[6,3])*sin(3*φ)) / ((r/R)**6) ) +  (  (legendre(4,6,sin0))*( (C[6,4])*cos(4*φ) + (S[6,4])*sin(4*φ)) / ((r/R)**6) )+  (  (legendre(5,6,sin0))*( (C[6,5])*cos(5*φ) + (S[6,5])*sin(5*φ)) / ((r/R)**6) ) +  (  (legendre(6,6,sin0))*( (C[6,6])*cos(6*φ) + (S[6,6])*sin(6*φ)) / ((r/R)**6) )
          + (J7 * legendre(0,7,sin0) / (r/R)**7) +    (  (legendre(1,7,sin0))*( (C[7,1])*cos(1*φ) + (S[7,1])*sin(1*φ)) / ((r/R)**7) ) + (  (legendre(2,7,sin0))*( (C[7,2])*cos(2*φ) + (S[7,2])*sin(2*φ)) / ((r/R)**7) ) + (  (legendre(3,7,sin0))*( (C[7,3])*cos(3*φ) + (S[7,3])*sin(3*φ)) / ((r/R)**7) ) +  (  (legendre(4,7,sin0))*( (C[7,4])*cos(4*φ) + (S[7,4])*sin(4*φ)) / ((r/R)**7) )+  (  (legendre(5,7,sin0))*( (C[7,5])*cos(5*φ) + (S[7,5])*sin(5*φ)) / ((r/R)**7) ) +  (  (legendre(6,7,sin0))*( (C[7,6])*cos(6*φ) + (S[7,6])*sin(6*φ)) / ((r/R)**7) ) +  (  (legendre(7,7,sin0))*( (C[7,7])*cos(7*φ) + (S[7,7])*sin(7*φ)) / ((r/R)**7) )
          + (J8 * legendre(0,8,sin0) / (r/R)**8)
          )


Fx = u.diff(x)
Fy = u.diff(y)
Fz = u.diff(z)
```

Finalmente el resultado obtenido es lo que se agrega al código original que fue entregado anteriormente en la entrega 6. La forma de implementar esto es muy sencilla ya que la funcion ***satelite_mejorado(z,t)*** utilizada en el odeint, queda reducida solamente a ensamblar el vector *"z punto"*. Otra mejora que se implementó fue definir las matrices de rotacion ***fuera de la función***, lo que aumentó el tiempo de resolución en un par de milisegundos. De esta manera la funcion queda como:

```
def satelite_mejorado(z,t):
    x = z[0:3]
    xp = z[3:6]
    zp = sp.zeros(6)
    zp[0:3] = xp
    zp[3:6] = Rtr @ (-Fg(x) - (((w*Rp2)@(xp)) + (((w**2)*Rpp)@(x))))
    return zp 
```

¿Pero donde queda la formula del potencial que mencione antes? Se puede ver que en la segunda mitad del vector zp, se incluye una funcion anidada, de manera que entregue justamente el vector de fuerzas necesitado en base al vector de posiciones ***"x"***. Este vector de fuerzas esta dado por la función ***Fg(z1)*** donde z1 representa el mismo vector de posiciones:


```
def Fg(z1):   
    x=z1[0]
    y=z1[1]
    z=z1[2]
    r2= z1[0]**2 + z1[1]**2 + z1[2]**2
    R = 6378136.3              # m         : Radio de la tierra
    μ = 398600441500000.000000 # m**3/s**2 : G*(masa de la tierra)
    J2 =   -0.10826360229840E-02  
    J3 =    0.25324353457544E-05   
    J4 =    0.16193312050719E-05  
    J5 =    0.22771610163688E-06  
    J6 =   -0.53964849049834E-06  
    J7 =    0.35136844210318E-06
    J8 =    0.20251871520885E-06
    
    Fx = x*μ*(J2*R**2*(3*z**2/(2*(r2)) - 0.5)/(r2) + J3*R**3*(5*z**3/(2*(r2)**(3/2)) -  {••••••••••••}  + 1.80773614603738e-6*y*cos(2*atan2(y, x))/(x**2 +
         y**2))/(r2))/sqrt(r2)
    Fy=y*μ*(J2*R**2*(3*z**2/(2*(r2)) - 0.5)/(r2) + J3*R**3*(5*z**3/(2*(r2)**(3/2)) - {••••••••••••} - 1.80773614603738e-6*x*cos(2*atan2(y, x))/(x**2 +
       y**2))/(r2))/sqrt(r2)
    Fz=z*μ*(J2*R**2*(3*z**2/(2*(r2)) - 0.5)/(r2) + J3*R**3*(5*z**3/(2*(r2)**(3/2)) - {••••••••••••} + 1.5745360427672e-6*cos(2*atan2(y, x)))/(r2))/sqrt(r2)
    
    return sp.array([Fx,Fy,Fz])
```

Esto se hizo así ya que si los términos eran incluidos dentro de la funcion  ***satelite_mejorado(z,t)*** , el tiempo de ejecucion era un par de milisegundos mayor. Se probaron los archivos 1 , 2 , 3 y 4 subidos EOF subidos en este repositorio (el archivo 1 corresponde al mío y los otros 3 son al azar) y se obtuvieron los siguientes tiempos y graficos de deriva respectivos:<br>
<br>

python3 entrega\ 7.py 1.EOF  21.79s user 0.57s system 29% cpu 1:15.65 total<br>
python3 entrega\ 7.py 2.EOF  21.09s user 0.32s system 22% cpu 1:34.97 total<br>
python3 entrega\ 7.py 3.EOF  22.82s user 0.59s system 78% cpu 29.717 total<br>
python3 entrega\ 7.py 4.EOF  23.21s user 0.58s system 88% cpu 26.975 total<br>
<br>



![alt text]( https://github.com/maxipoblete/MCOC2020-P1/blob/master/E7%201.png )

![alt text](  https://github.com/maxipoblete/MCOC2020-P1/blob/master/E7%202.png  )

![alt text](  https://github.com/maxipoblete/MCOC2020-P1/blob/master/E7%203.png  )

![alt text](  https://github.com/maxipoblete/MCOC2020-P1/blob/master/E7%204.png  )
<br>
<br>
Se puede apreciar que la mejora no es tanta la mejora en la deriva y al final el tiempo de ejecucion es mucho mayor....Pero esto debe verse al final cuando en el concurso se corran multiples archivos EOF, ya que cada codigo es distinto y cada archivo es distino. Hay que aumentar la muestra para obtener una vision mas general.



