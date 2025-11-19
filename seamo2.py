# -*- coding: utf-8 -*-
"""
seamo2.py
Implementación del algoritmo SEAMO2 para el URTP,
siguiendo el paper CEC 2013 (Mumford) y sus Fig. 3, 4 y 5.
Evaluación del coste del pasajero (CP) acelerada con Floyd–Warshall
(una única ejecución APSP por individuo).
"""

import random #para toda la aleatoriedad del algoritmo
import heapq  #no lo uso al final, pues para calcular el CP uso F-W en vez de Dijkstra
from dataclasses import dataclass, field #para definir contenedores inmutables por diseño y con menos código
from typing import List, Dict, Set, Tuple, Optional
import collections

INF = 1e12 # Un “infinito” numérico grande para inicializar distancias en Floyd–Warshall y detectar pares inalcanzables

# -------------------------------------------------------------
# 1) Estructuras de red y rutas
# -------------------------------------------------------------

@dataclass
class GrafoTransporte: #GrafoTransporte = red de transporte urbano
    n: int #n=numero de nodos de la red
    ady: Dict[int, Dict[int, float]] = field(default_factory=dict)
    demanda: List[List[float]] = field(default_factory=list)
    penalizacion_transbordo: float = 5.0

    def __post_init__(self):
        if not self.ady:
            self.ady = {i: {} for i in range(self.n)}
        if not self.demanda:
            self.demanda = [[0.0]*self.n for _ in range(self.n)]

    def agregar_arista(self, u: int, v: int, w: float):
        # grafo no dirigido y tiempos simétricos
        self.ady[u][v] = float(w)
        self.ady[v][u] = float(w)

@dataclass
class Ruta:
    nodos: List[int]

    def aristas(self):
        return [(self.nodos[i], self.nodos[i+1]) for i in range(len(self.nodos)-1)]

# --- Utilidades para evitar duplicados de soluciones ---
def canonical_route(route):
    """Devuelve la orientación canónica de una ruta: la menor lexicográficamente entre r y reversed(r)."""
    r = tuple(route)
    rr = tuple(reversed(route))
    return r if r < rr else rr

def canonical_signature(rutas):
    """
    Firma canónica de un conjunto de rutas:
    - normaliza la orientación de cada ruta (independiente de invertirla)
    - ordena las rutas (independiente del orden en el conjunto)
    """
    norm = [canonical_route(r.nodos) if hasattr(r, "nodos") else canonical_route(r) for r in rutas]
    norm.sort()
    return tuple(norm)


@dataclass
class ConjuntoRutas:
    rutas: List[Ruta]

    def nodos_usados(self) -> Set[int]:
        usados = set()
        for r in self.rutas:
            usados.update(r.nodos)
        return usados

    def copy(self) -> "ConjuntoRutas":
        return ConjuntoRutas([Ruta(r.nodos[:]) for r in self.rutas])

# -------------------------------------------------------------
# 2) Red de tránsito y evaluación (CP con Floyd–Warshall) + CO
# -------------------------------------------------------------

class RedTransito:
    """
    Cada ocurrencia de un nodo (aparición en alguna ruta) es un vértice distinto.
    - Aristas de transporte: entre ocurrencias consecutivas en una ruta, con tiempo de viaje.
    - Aristas de transbordo: conectan ocurrencias del mismo nodo entre rutas, con penalización fija.
    """
    def __init__(self, G: GrafoTransporte, R: ConjuntoRutas):
        self.G = G
        self.R = R
        self.vertice_de = {}  # (indice_ruta, pos_en_ruta) -> id_vertice
        self.ocurrencias = collections.defaultdict(list)  # nodo -> [ids_ocurrencia]
        self.ady = {}  # lista de adyacencia en red de tránsito: u -> [(v,w),...]
        self.n_vertices = 0
        self._construir()

    def _construir(self):
        # Crear vértices: una ocurrencia por posición en cada ruta
        for a, ruta in enumerate(self.R.rutas):
            for k, nodo in enumerate(ruta.nodos):
                vid = self.n_vertices
                self.vertice_de[(a, k)] = vid
                self.ocurrencias[nodo].append(vid)
                self.n_vertices += 1

        # Adjacencias
        self.ady = {i: [] for i in range(self.n_vertices)}

        # Aristas de transporte (ambas direcciones: ida/vuelta en la línea)
        for a, ruta in enumerate(self.R.rutas):
            for k in range(len(ruta.nodos) - 1):
                u_node = ruta.nodos[k]
                v_node = ruta.nodos[k + 1]
                w = self.G.ady[u_node][v_node]
                u = self.vertice_de[(a, k)]
                v = self.vertice_de[(a, k + 1)]
                self.ady[u].append((v, w))
                self.ady[v].append((u, w))

        # Aristas de transbordo entre ocurrencias del mismo nodo
        p = self.G.penalizacion_transbordo
        for _, vids in self.ocurrencias.items():
            L = len(vids)
            if L >= 2:
                for i in range(L):
                    for j in range(i+1, L):
                        u, v = vids[i], vids[j]
                        self.ady[u].append((v, p))
                        self.ady[v].append((u, p))

    def floyd_warshall(self):
        """APSP (todos los pares) sobre la red de tránsito: una sola vez por individuo."""
        n = self.n_vertices
        dist = [[INF]*n for _ in range(n)]
        for u in range(n):
            dist[u][u] = 0.0
            for (v, w) in self.ady[u]:
                if w < dist[u][v]:
                    dist[u][v] = w
        # FW triple bucle
        for k in range(n):
            dk = dist[k]
            for i in range(n):
                dik = dist[i][k]
                if dik == INF:
                    continue
                di = dist[i]
                base = dik
                for j in range(n):
                    nd = base + dk[j]
                    if nd < di[j]:
                        di[j] = nd
        return dist

def coste_pasajero(G: GrafoTransporte, R: ConjuntoRutas) -> float:
    # PRECHECK: si algún nodo con demanda no aparece en rutas, inviable (evita correr FW inútilmente)
    T = RedTransito(G, R)
    for i in range(G.n):
        for j in range(G.n):
            dij = G.demanda[i][j]
            if i != j and dij > 0:
                if not T.ocurrencias.get(i) or not T.ocurrencias.get(j):
                    return INF

    # Floyd–Warshall una sola vez (como manda el paper)
    dist = T.floyd_warshall()

    num, den = 0.0, 0.0
    for i in range(G.n):
        fuentes = T.ocurrencias.get(i, [])
        if not fuentes:
            return INF
        for j in range(G.n):
            dij = G.demanda[i][j]
            if i == j or dij <= 0:
                continue
            destinos = T.ocurrencias.get(j, [])
            if not destinos:
                return INF

            mejor = INF
            for s in fuentes:
                ds = dist[s]
                for t in destinos:
                    if ds[t] < mejor:
                        mejor = ds[t]
            if mejor >= INF/2:
                return INF  # O–D inalcanzable → inviable

            num += dij * mejor
            den += dij
    return num/den if den > 0 else INF




def coste_operador(G: GrafoTransporte, R: ConjuntoRutas) -> float:
    total = 0.0
    for r in R.rutas:
        for (u, v) in r.aristas():
            total += G.ady[u][v]
    return total

# -------------------------------------------------------------
# 3) Inicialización y reparación (Fig. 4 y Fig. 5)
# -------------------------------------------------------------

@dataclass
class Restricciones:
    r: int
    MIN: int
    MAX: int

def generar_rutaset(G: GrafoTransporte, C: Restricciones) -> Optional[ConjuntoRutas]:
    """
    Fig. 4 (Mumford): Inicialización por crecimiento de rutas nodo-a-nodo,
    evitando ciclos, usando un 'reverse' si se atasca (un intento extra), y reparación posterior.
    """
    elegidos: Set[int] = set()
    rutas: List[Ruta] = []

    for cuenta in range(1, C.r + 1):
        L = random.randint(C.MIN, C.MAX)

        # 1) Semilla:
        #    - 1ª ruta (o sin elegidos): nodo con grado > 0 si es posible
        #    - siguientes: preferir un nodo ya usado con grado > 0
        if cuenta == 1 or not elegidos:
            cand0 = [u for u, vec in G.ady.items() if len(vec) > 0]
            i = random.choice(cand0) if cand0 else random.randrange(G.n)
        else:
            candE = [u for u in elegidos if len(G.ady[u]) > 0]
            i = random.choice(candE) if candE else random.randrange(G.n)

        ruta = [i]
        elegidos.add(i)

        # 2) Crecimiento con reverse (un intento desde el otro extremo)
        invertido_una_vez = False
        while len(ruta) < L:
            u = ruta[-1]
            candidatos = [v for v in G.ady[u] if v not in ruta]  # evita ciclos internos
            if candidatos:
                v = random.choice(candidatos)
                ruta.append(v)
                elegidos.add(v)
            else:
                # si la ruta tiene >1 nodo, invertir y probar una sola vez desde el otro extremo
                if not invertido_una_vez and len(ruta) > 1:
                    ruta.reverse()
                    invertido_una_vez = True
                    continue  # reintentar crecimiento desde el nuevo extremo
                break  # atascado incluso tras el reverse o ruta de tamaño 1

        rutas.append(Ruta(ruta))

    rutaset = ConjuntoRutas(rutas)

    # 3) Reparación de cobertura (Fig. 5)
    if len(rutaset.nodos_usados()) < G.n:
        rutaset = reparar_rutaset(G, rutaset, C)
        if rutaset is None:
            return None

    return rutaset


def reparar_rutaset(G: GrafoTransporte, R: ConjuntoRutas, C: Restricciones) -> Optional[ConjuntoRutas]:
    """
    Fig. 5 (Mumford): Inserta nodos faltantes en extremos de rutas si están
    adyacentes y no superan MAX. Si no puede colocar alguno, devuelve None.
    """
    faltan = [v for v in range(G.n) if v not in R.nodos_usados()]
    rutas = [Ruta(r.nodos[:]) for r in R.rutas]
    for v in faltan:
        colocado = False
        indices = list(range(len(rutas)))
        random.shuffle(indices)
        for idx in indices:
            if len(rutas[idx].nodos) >= C.MAX:
                continue
            izq = rutas[idx].nodos[0]
            der = rutas[idx].nodos[-1]
            if v in G.ady[izq]:
                rutas[idx].nodos = [v] + rutas[idx].nodos
                colocado = True
                break
            if v in G.ady[der]:
                rutas[idx].nodos.append(v)
                colocado = True
                break
        if not colocado:
            return None
    return ConjuntoRutas(rutas)

# -------------------------------------------------------------
# 4) Operadores: cruce y mutación (según espíritu del paper)
# -------------------------------------------------------------

def cruce(G: GrafoTransporte, P1: ConjuntoRutas, P2: ConjuntoRutas, C: Restricciones) -> ConjuntoRutas:
    import random

    def elegir_mejor(eligibles, usados):
        best, best_score = None, -1.0
        for r in eligibles:
            ruta = r.nodos
            nuevos = sum(1 for x in ruta if x not in usados)
            score = nuevos / max(1, len(ruta))
            if score > best_score:
                best, best_score = r, score
        return best

    rutas1 = P1.copy().rutas[:]; random.shuffle(rutas1)
    rutas2 = P2.copy().rutas[:]; random.shuffle(rutas2)
    hijo, usados = [], set()

    # Semilla de P1
    if not rutas1: rutas1, rutas2 = rutas2, rutas1
    semilla = rutas1.pop()
    hijo.append(Ruta(semilla.nodos[:])); usados.update(semilla.nodos)

    # Alternar padres, eligiendo SOLO rutas que compartan al menos un nodo con lo ya elegido
    turn = 2
    while len(hijo) < C.r:
        candidatos = rutas2 if turn == 2 else rutas1
        elegibles = [r for r in candidatos if any(x in usados for x in r.nodos)]
        if elegibles:
            elegido = elegir_mejor(elegibles, usados)
            candidatos.remove(elegido)
            hijo.append(Ruta(elegido.nodos[:])); usados.update(elegido.nodos)
            turn = 1 if turn == 2 else 2
        else:
            other = rutas1 if turn == 2 else rutas2
            elegibles_other = [r for r in other if any(x in usados for x in r.nodos)]
            if elegibles_other:
                elegido = elegir_mejor(elegibles_other, usados)
                other.remove(elegido)
                hijo.append(Ruta(elegido.nodos[:])); usados.update(elegido.nodos)
            else:
                pool = candidatos if candidatos else other
                if not pool: break
                elegido = elegir_mejor(pool, usados)
                pool.remove(elegido)
                hijo.append(Ruta(elegido.nodos[:])); usados.update(elegido.nodos)

    # Relleno extremo si faltan rutas
    while len(hijo) < C.r:
        u = random.randrange(G.n)
        vecinos = [v for v in G.ady[u]]
        if vecinos:
            v = random.choice(vecinos); hijo.append(Ruta([u, v])); usados.update([u, v])
        else:
            hijo.append(Ruta([u]))
    return ConjuntoRutas(hijo[:C.r])






def mutacion(G: GrafoTransporte, R: ConjuntoRutas, C: Restricciones, p_add_del: float = 0.5) -> ConjuntoRutas:
    """
    Mutación tipo add/del de extremos (evitar “bloating”), en línea con Mumford.
    """
    rutas = [Ruta(r.nodos[:]) for r in R.rutas]
    I = random.randint(1, max(1, (C.r * C.MAX) // 2))

    if random.random() < p_add_del:
        # ADD: añadir en extremos si hay adyacencia y no hay ciclo
        intentos = 0
        while intentos < I:
            idx = random.randrange(len(rutas))
            rt = rutas[idx]
            if len(rt.nodos) >= C.MAX:
                intentos += 1
                continue
            extremos = [rt.nodos[0], rt.nodos[-1]]
            end = random.choice(extremos)
            candidatos = [v for v in G.ady[end] if v not in rt.nodos]
            if not candidatos:
                intentos += 1
                continue
            v = random.choice(candidatos)
            if end == rt.nodos[0]:
                rutas[idx] = Ruta([v] + rt.nodos)
            else:
                rutas[idx] = Ruta(rt.nodos + [v])
            intentos += 1
    else:
        # DEL: quitar de extremos si el nodo aparece en otras rutas (para mantener cobertura)
        intentos = 0
        while intentos < I:
            idx = random.randrange(len(rutas))
            rt = rutas[idx]
            if len(rt.nodos) <= C.MIN:
                intentos += 1
                continue
            # contar coberturas globales
            cnt = collections.Counter()
            for r in rutas:
                cnt.update(r.nodos)
            opciones = []
            for endpos in [0, -1]:
                node = rt.nodos[endpos]
                if cnt[node] >= 2:  # aparece en otro sitio
                    opciones.append(endpos)
            if not opciones:
                intentos += 1
                continue
            endpos = random.choice(opciones)
            if endpos == 0:
                rutas[idx] = Ruta(rt.nodos[1:])
            else:
                rutas[idx] = Ruta(rt.nodos[:-1])
            intentos += 1

    return ConjuntoRutas(rutas)

# -------------------------------------------------------------
# 5) Algoritmo SEAMO2 (Fig. 3)
# -------------------------------------------------------------

@dataclass
class Individuo:
    rutaset: ConjuntoRutas
    cp: float
    co: float

    def domina(self, otro: "Individuo") -> bool:
        return (self.cp <= otro.cp and self.co <= otro.co) and (self.cp < otro.cp or self.co < otro.co)

def crear_individuo(G: GrafoTransporte, R: ConjuntoRutas) -> Individuo:
    return Individuo(R, coste_pasajero(G, R), coste_operador(G, R))

def seamo2(G: GrafoTransporte, C: Restricciones, tam_pob: int = 50, generaciones: int = 200, seed: Optional[int] = None) -> List[Individuo]:
    if seed is not None:
        random.seed(seed)

    # Población inicial
    poblacion: List[Individuo] = []
    intentos = 0
    while len(poblacion) < tam_pob and intentos < tam_pob * 50:
        R0 = generar_rutaset(G, C)
        intentos += 1
        if R0 is None:
            continue
        poblacion.append(crear_individuo(G, R0))
    if not poblacion:
        raise RuntimeError("No se pudo crear población inicial. Revisa restricciones y grafo.")

    # Mejores por objetivo (seguimiento)
    mejor_cp = min(poblacion, key=lambda x: x.cp)
    mejor_co = min(poblacion, key=lambda x: x.co)

    # Bucle evolutivo
    for _ in range(generaciones):
        for i in range(len(poblacion)):
            p1 = poblacion[i]
            j = random.randrange(len(poblacion)-1)
            if j >= i:
                j += 1
            p2 = poblacion[j]

            # Crossover → reparación → mutación
            hijo_R = cruce(G, p1.rutaset, p2.rutaset, C)
            hijo_R = reparar_rutaset(G, hijo_R, C) or hijo_R
            hijo_R = mutacion(G, hijo_R, C)
            # --- Deduplicar rutas del hijo (usando tu canonical_route) y rellenar sin duplicar ---
            unicas = []
            vistos = set()

            for r in hijo_R.rutas:
                sig = canonical_route(r.nodos)  # trata [1,2,3] y [3,2,1] como la misma ruta
                if sig not in vistos:
                    unicas.append(Ruta(list(sig)))  # normalizamos orientación
                    vistos.add(sig)

            # Relleno: genera rutas cortas nuevas que no estén ya en 'vistos'
            guard = 0
            while len(unicas) < C.r and guard < 10_000:
                # ancla en nodos ya presentes si es posible (mejora conectividad), si no al azar
                if unicas:
                    ancla = {x for rr in unicas for x in rr.nodos}
                    u = random.choice(tuple(ancla))
                else:
                    u = random.randrange(G.n)

                vecinos = list(G.ady[u].keys())
                if not vecinos:
                    guard += 1
                    continue

                v = random.choice(vecinos)
                cand = [u, v]

                # extender hasta MIN evitando ciclos cortos
                while len(cand) < C.MIN:
                    ult = cand[-1]
                    opciones = [w for w in G.ady[ult].keys() if w not in cand]
                    if not opciones:
                        break
                    cand.append(random.choice(opciones))

                if C.MIN <= len(cand) <= C.MAX:
                    sig = canonical_route(cand)
                    if sig not in vistos:
                        unicas.append(Ruta(list(sig)))
                        vistos.add(sig)

                guard += 1

            hijo_R = ConjuntoRutas(unicas)
            # -------------------------------------------------------------------------------

            # Reparación final básica: cobertura
            if len(hijo_R.rutas) != C.r:
                continue
            if len(hijo_R.nodos_usados()) < G.n:
                fix = reparar_rutaset(G, hijo_R, C)
                if fix is None:
                    continue
                hijo_R = fix

            hijo = crear_individuo(G, hijo_R)

            # Reglas de reemplazo estilo SEAMO2 (simplificadas pero fieles al espíritu)
            if hijo.domina(p1):
                poblacion[i] = hijo
            elif hijo.domina(p2):
                poblacion[j] = hijo
            else:
                # Insertar si mejora algún extremo del frente
                if hijo.cp < mejor_cp.cp or hijo.co < mejor_co.co:
                    if (p1.cp + p1.co) > (p2.cp + p2.co):
                        poblacion[i] = hijo
                    else:
                        poblacion[j] = hijo

            # Actualizar marcadores
            if poblacion[i].cp < mejor_cp.cp:
                mejor_cp = poblacion[i]
            if poblacion[j].co < mejor_co.co:
                mejor_co = poblacion[j]

    # Extraer no dominados finales
    # Extraer no dominados finales (únicos por firma canónica)
    frente = []
    vistas = set()
    for ind in poblacion:
        if any(otro.domina(ind) for otro in poblacion if otro is not ind):
            continue
        # Normaliza orientación de cada ruta para una salida estable
        for r in ind.rutaset.rutas:
            r.nodos = list(canonical_route(r.nodos))
        sig = canonical_signature(ind.rutaset.rutas)
        if sig in vistas:
            continue
        vistas.add(sig)
        frente.append(ind)
    return frente

def indicadores_demanda(G, rutaset):
    """
    Calcula d0, d1, d2, d_un (en %) para un conjunto de rutas.
    Modelo: transbordos mínimos entre rutas, con máximo permitido de 2.
    Si requiere >2 transbordos o es inalcanzable, cuenta como demanda insatisfecha.
    """
    rutas = [list(r.nodos) if hasattr(r, "nodos") else list(r) for r in rutaset.rutas]

    # índice: nodo -> set(indices de rutas que lo contienen)
    rutas_por_nodo = [set() for _ in range(G.n)]
    for i, r in enumerate(rutas):
        for u in r:
            rutas_por_nodo[u].add(i)

    # grafo de rutas: arista si comparten al menos un nodo
    adj_rutas = [set() for _ in range(len(rutas))]
    for u in range(G.n):
        lst = list(rutas_por_nodo[u])
        for a in range(len(lst)):
            for b in range(a + 1, len(lst)):
                ra, rb = lst[a], lst[b]
                adj_rutas[ra].add(rb)
                adj_rutas[rb].add(ra)

    tot_dem = 0.0
    dem0 = dem1 = dem2 = 0.0

    for o in range(G.n):
        for d in range(G.n):
            dem = G.demanda[o][d]
            if dem <= 0 or o == d:
                continue
            tot_dem += dem

            rutas_o = rutas_por_nodo[o]
            rutas_d = rutas_por_nodo[d]

            if not rutas_o or not rutas_d:
                continue  # no cubierto -> cuenta en d_un

            if rutas_o & rutas_d:
                dem0 += dem  # mismo recorrido: 0 transbordos
                continue

            # 1 transbordo: existe ruta ro y rd que se crucen
            un_transfer = False
            for ro in rutas_o:
                if any(rd in adj_rutas[ro] for rd in rutas_d):
                    dem1 += dem
                    un_transfer = True
                    break
            if un_transfer:
                continue

            # 2 transbordos: ro -- r? -- rd
            dos_transfer = False
            for ro in rutas_o:
                for rmid in adj_rutas[ro]:
                    if any(rd in adj_rutas[rmid] for rd in rutas_d):
                        dem2 += dem
                        dos_transfer = True
                        break
                if dos_transfer:
                    break
            # si no hay 0/1/2 -> queda como insatisfecha

    if tot_dem == 0:
        return 0.0, 0.0, 0.0, 100.0

    d0 = 100.0 * dem0 / tot_dem
    d1 = 100.0 * dem1 / tot_dem
    d2 = 100.0 * dem2 / tot_dem
    d_un = max(0.0, 100.0 - (d0 + d1 + d2))
    return d0, d1, d2, d_un


# -------------------------------------------------------------
# 6) Demo mínimo (red toy) — opcional para comprobar instalación
# -------------------------------------------------------------

def red_ejemplo() -> GrafoTransporte:
    G = GrafoTransporte(6)
    aristas = [(0,1,2),(1,2,2),(2,3,2),(3,4,2),(4,5,2),(5,0,2),(1,4,3)]
    for u,v,w in aristas:
        G.agregar_arista(u,v,w)
    for i in range(6):
        for j in range(6):
            if i != j:
                G.demanda[i][j] = 1.0
    return G

if __name__ == "__main__":
    # Prueba rápida con red toy
    random.seed(1)
    G = red_ejemplo()
    C = Restricciones(r=3, MIN=3, MAX=5)
    frente = seamo2(G, C, tam_pob=10, generaciones=30)
    print("Soluciones no dominadas (toy):")
    for ind in frente:
        if ind.cp < 1e9:
            rutas = [[x for x in r.nodos] for r in ind.rutaset.rutas]
            print(f"CP={ind.cp:.2f}, CO={ind.co:.2f}, rutas={rutas}")
