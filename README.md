# Turbulencias

Proyecto Metodos de Simulación

### Por favor interntar comentar lo maximo posible el codigo que siban

## Carpetas

En este proyecto, se encuentran las siguientes carpetas:

1. `src`: Esta carpeta contiene el código fuente del proyecto, donde se encuentran las diferentes etapas de la simulacion, cada una con un codigo, tambien estan los codigos para anlizar los datos.

2. `uotp`: En esta carpeta se almacenan los resultados de las diferentes simulaciones, tanto graficas, gif como datos crudos

3. `run`: Carpeta dedicada unicamente a correr los programas de forma mas eficiente

## Estructura del Proyecto

El proyecto está dividido en las siguientes etapas:

1. _Movimiento Browniano del Disco_
2. _Disco con Velocidad Angular y Fuerza de Magnus_
3. _Disco con Fuerza Central_
4. _Simulación con Múltiples Discos_
5. _Simulación de Fluidos con Lattice Boltzmann (LB)_
6. _Analisis de datos_

### 1. Movimiento Browniano del Disco

En esta etapa, se simula el movimiento browniano de un disco. El código de esta simulación es en src/Disco_difusion.cpp

- _Input_: El archivo de configuración acepta un parámetro que puede ser 0 o 1.

  - 0: Imprime los datos crudos.
  - 1: Imprime los datos en formato compatible con gnuplot para generar un gif.

- _Output_:
  - Datos crudos de la simulación o archivos preparados para gnuplot, dependiendo del input especificado.

### 2. Disco con Velocidad Angular y Fuerza de Magnus

En esta etapa, se agrega una velocidad angular y la fuerza de Magnus al disco. El código de esta simulación es en src/Disco_Magnus.cpp

- _Input_: El archivo de configuración acepta un parámetro que puede ser 0 o 1.

  - 0: Imprime los datos crudos.
  - 1: Imprime los datos en formato compatible con gnuplot para generar un gif.

- _Output_:
  - Datos crudos de la simulación o archivos preparados para gnuplot, dependiendo del input especificado.

### 3. Disco con Fuerza Central

En esta etapa, se añade una fuerza central que actúa sobre el disco. El código de esta simulación es en src/Disco_Fcentral.cpp

- _Input_: El archivo de configuración acepta un parámetro que puede ser 0 o 1.

  - 0: Imprime los datos crudos.
  - 1: Imprime los datos en formato compatible con gnuplot para generar un gif.

- _Output_:
  - Datos crudos de la simulación o archivos preparados para gnuplot, dependiendo del input especificado.

### 4. Simulación con Múltiples Discos

En esta etapa, se simula el comportamiento de múltiples discos interactuando entre sí. El código de esta por decidirse.

- _Input_: El archivo de configuración acepta un parámetro que puede ser 0 o 1.

  - 0: Imprime los datos crudos.
  - 1: Imprime los datos en formato compatible con gnuplot para generar un gif.

- _Output_:
  - Datos crudos de la simulación o archivos preparados para gnuplot, dependiendo del input especificado.

### 5. Simulación de Fluidos con Lattice Boltzmann (LB)

En esta etapa final, se utiliza el método de Lattice Boltzmann para simular el comportamiento de fluidos. El código de esta por decidirse.

- _Input_:

  - Por ver

- _Output_:
  - Por ver

### 6. Analisis de datos

Existen dos codigos para esta etapa, src/Analisis.ipynb es un notebook en donde se ejecuta todo, incluyendo algunas pruebas. src/Paper.py idelamente solo debe ejecutar lo importante

- _Input_:

  - datos crudos

- _Output_:
  - Graficas

## Cómo Ejecutar las Simulaciones

1. Idealmente las simulaciones se corren desde los scripts de la shel .sh
2. Se pueden correr los analisis de datos de dos formas difernetes, ya sea celda por celda en el notebook, o en el script de python para ejecutar todo de una (cuando este)
