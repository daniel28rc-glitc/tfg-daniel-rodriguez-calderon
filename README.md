# TFG - Daniel Rodríguez Calderón

Este repositorio contiene los scripts implementados en MATLAB para el Trabajo de Fin de Grado de Daniel Rodríguez Calderón. Los códigos incluyen diversas simulaciones y cálculos enfocados en el sistema de ecuaciones diferenciales estocásticas, la ecuación en derivadas parciales asociada o las soluciones semi-analíticas del modelo de Heston.

## Archivos del repositorio

- **movbrown.m**: Script enfocado en la simulación matemática y generación de múltiples trayectorias de movimientos brownianos.
- **mallaEDP.m**: Código que dibuja de forma generalizada el dominio truncado de evalución del método EDP.
- **EDP.m**: Archivo que desarrolla la valoración de opciones europeas mediante el método EDP.
- **EDE.m**: Script que implementa la valoración de opciones europeas mediante el método EDE.
- **SA.m**: Código que realiza la valoración de opciones europeas mediante el método SA.
- **ordEDP.m**: Archivo que evalúa los órdenes de convergencia espacio-temporales del método EDP.
- **ordEDE.m**: Script que calcula los órdenes de convergencia fuertes y débiles del método EDE.
- **compLp.m**: Código que compara los diferentes métodos de valoración: EDP, EDE y SA.

## Requisitos de uso

Para ejecutar correctamente los archivos de este repositorio, es necesario:
1. Contar con el software MATLAB instalado en su equipo.
2. Descargar los archivos `.m` en un mismo directorio.
3. Ejecutarlos desde ese mismo entorno de MATLAB.
