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


