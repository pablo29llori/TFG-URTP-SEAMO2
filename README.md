# 🚍 Trabajo Fin de Grado – El problema de diseño de rutas en el tránsito urbano.  
## Grado en Matemáticas – Universidad de Oviedo
## Autor: Pablo Llorián González
## Defensa: 18 de noviembre de 2025


Este repositorio recoge el material principal de mi Trabajo Fin de Grado, titulado **“El problema de diseño de rutas en el tránsito urbano”**, en el que se estudia el **Urban Transit Routing Problem (URTP)** desde una perspectiva de optimización combinatoria y análisis multiobjetivo.  

El objetivo es diseñar un conjunto de rutas de transporte urbano que dé servicio a una red de nodos (paradas), equilibrando el compromiso entre:  

- **Coste del operador (CO)**: longitud total de las rutas.  
- **Coste del pasajero (CP)**: tiempo medio de viaje, incluyendo transbordos y penalizaciones asociadas.  

Para ello se sigue el enfoque del paper de Christine Mumford (CEC 2013), adaptando e implementando el algoritmo evolutivo **SEAMO2** y aplicándolo sobre diferentes redes de transporte de prueba.  

---

## 📁 Contenido del repositorio

### **1. Memoria del TFG**  
Documento completo del Trabajo Fin de Grado, donde se desarrollan:  

- La formulación del **URTP** como problema de optimización en redes.  
- La base teórica de y algoritmos para optimización.  
- La descripción de los objetivos **CP** (coste del pasajero) y **CO** (coste del operador).  
- El resumen del enfoque propuesto en el artículo de Mumford (CEC 2013).  
- El diseño detallado del algoritmo **SEAMO2** adaptado al URTP.  
- La aceleración del cálculo del CP mediante una única ejecución de **Floyd–Warshall** por individuo.  
- La descripción de los conjuntos de datos utilizados (red AVE y red Mumford2).  
- Los experimentos numéricos, análisis de frentes de Pareto y métricas de calidad (incluyendo indicadores de demanda servida con 0, 1 o 2 transbordos).  
- Discusión de resultados y posibles líneas futuras de trabajo.  

📄 `TFG-Diseno-Rutas-Transito-Urbano.pdf`.  

---

### **2. Implementación del algoritmo SEAMO2 en Python**  

El archivo `seamo2.py` contiene una implementación completa del algoritmo **SEAMO2** para el URTP, siguiendo la estructura y el espíritu del artículo original:  

- Definición de la red de transporte urbano como grafo ponderado.  
- Representación de soluciones como conjuntos de rutas que cubren todos los nodos bajo restricciones de longitud y número de rutas.  
- Construcción de la **red de tránsito** (apariciones de nodos en rutas, aristas de viaje y transbordo).  
- Cálculo del **coste del pasajero (CP)** mediante **Floyd–Warshall** (APSP) en la red de tránsito.  
- Cálculo del **coste del operador (CO)** como suma de las longitudes de las rutas.  
- Generación de soluciones iniciales y reparación de cobertura.  
- Operadores evolutivos:  
  - Cruce de conjuntos de rutas.  
  - Mutación mediante ampliación/recorte de extremos.  
  - Reparaciones adicionales para mantener cobertura y restricciones.  
- Algoritmo principal **SEAMO2**, con:  
  - Comparación por dominancia (CP, CO).  
  - Mantenimiento de una población de soluciones no dominadas.  
  - Extracción del frente de Pareto final.  

El script incluye además una **red de ejemplo pequeña** para comprobar que la implementación funciona correctamente.  

📄 `seamo2.py`.  

---

### **3. Conjuntos de datos**  

En la carpeta `datos/` se incluyen dos redes de prueba, cada una descrita mediante tres ficheros `.csv` (nodos, arcos y demanda):  

📁 `datos/`  

- **Red AVE.**  
  - `ave_nodes.csv`: nodos de la red (identificadores y, en su caso, información asociada).  
  - `ave_links.csv`: arcos de la red con sus longitudes/tiempos.  
  - `ave_demand.csv`: matriz (o lista) de demanda entre pares de nodos.  

- **Red Mumford2.**  
  - `mumford2_nodes.csv`.  
  - `mumford2_links.csv`.  
  - `mumford2_demand.csv`.  

Estos ficheros permiten reconstruir la red de transporte urbano y la demanda asociada para reproducir los experimentos realizados en el TFG.  

---

## 👤 Autor y contribución personal 

Soy **Pablo Llorián González**, autor del Trabajo Fin de Grado y de la implementación en Python del algoritmo **SEAMO2** incluida en este repositorio. Mis principales contribuciones son:  

- Estudio teórico y formulación del **URTP** y sus objetivos.  
- Adaptación del algoritmo **SEAMO2** al contexto del problema de diseño de rutas en tránsito urbano.  
- Implementación en Python del modelo de red, evaluación de CP y CO, generación de soluciones y operadores evolutivos.  
- Preparación y tratamiento de los conjuntos de datos utilizados en los experimentos.  
- Diseño, ejecución y análisis de los experimentos numéricos.  
- Redacción completa de la memoria del TFG.  

---

## 🛠️ Tecnologías utilizadas

- **Python** para la implementación del algoritmo SEAMO2.  
- **Floyd–Warshall** para el cálculo de caminos mínimos en la red de tránsito.  
- **Pandas / CSV** para el tratamiento de datos de nodos, enlaces y demanda.  
- **LaTeX** para la redacción de la memoria del TFG.  
- **GitHub** como repositorio y portfolio del proyecto.  

---

## 📅 Fecha

DepósitO: 7 de noviembre de 2025.
Defensa: 18 de noviembre de 2025.  
