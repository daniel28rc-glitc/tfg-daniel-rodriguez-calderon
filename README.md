# TFG - Daniel Rodríguez Calderón

Este repositorio contiene los scripts implementados en MATLAB para el Trabajo de Fin de Grado de Daniel Rodríguez Calderón. Los códigos incluyen diversas simulaciones y cálculos enfocados en ecuaciones diferenciales estocásticas, derivadas parciales o soluciones semi-analíticas.

## Archivos del repositorio

- **MovBrowniano.m**: Script enfocado en la simulación matemática y generación de múltiples trayectorias de movimientos brownianos, permitiendo visualizar distintas trayectorias.
- **EDE.m**: Implementa la valoración de opciones europeas en el Modelo de Heston utilizando el método de Euler-Maruyama combinado con simulaciones de Monte Carlo.
- **EDP.m**: Desarrolla la valoración de opciones europeas en el Modelo de Heston mediante un método en Ecuaciones en Derivadas Parciales (EDP), utilizando Euler explícito con diferencias finitas.
- **SA.m**: Realiza la valoración de opciones europeas en el Modelo de Heston utilizando un método semi-analítico basado en integración numérica.
- **compLp.m**: Compara diferentes métodos de valoración en el Modelo de Heston mediante el cálculo de la norma $L^p(Q)$ utilizando el teorema de Feynman-Kac.
- **ordEDE.m**: Script auxiliar que calcula los órdenes de convergencia del método de Euler-Maruyama junto con Monte Carlo.
- **ordEDP.m**: Script auxiliar que evalúa los órdenes de convergencia del método en EDP mediante diferencias finitas.
- **mallaEDP.m**: Dibuja de forma general, el esquema discreto de evalución del método en EDP.
## Requisitos de uso

Para ejecutar correctamente los archivos de este repositorio, es necesario:
1. Contar con el software MATLAB instalado en su equipo.
2. Descargar los archivos `.m` en un mismo directorio.
3. Ejecutarlos desde ese mismo entorno de MATLAB.
