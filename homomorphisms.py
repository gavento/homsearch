"""
Tools for finding graph homomorphisms and cores.

Loads a Cython module for the actual search routines.
This module itself contains mostly auxiliary tools:
* Conversion Graph <-> raw CGraph
* Heuristics for branching order
* Finding and verifying homomorphism cores
"""

import logging
from random import shuffle

from sage.graphs.graph import Graph
from sage.misc.preparser import load
from sage.graphs.base.dense_graph import DenseGraph


log = logging.getLogger("hom")
log.setLevel(logging.INFO)

# Initialize logging if not already initialized
if not logging.root.handlers:
  logging.basicConfig(format="%(asctime)s %(levelname)s [%(name)s]: %(msg)s", datefmt="%Y-%m-%d %H:%M:%S")

try:
  load("homomorphisms_c.pyx", locals())
except Exception as e:
  log.fatal("homomorphisms_c.pyx failed to load: %s", e)


### Heurisrics for branching vertex order


def degree_within(G, v, w):
  "Return the number of neighbors of `v` in `w`."
  return len(set(G.neighbors(v)).intersection(set(w)))


def sort_vertices_by_degree(G, vs, within=None):
  """
  Stable inplace sort of vs by degree in G, optionally
  only considering neighbors in `within`.
  """
  if within is None:
    vs.sort(cmp=lambda u,v: cmp(G.degree(u), G.degree(v)))
  else:
    vs.sort(cmp=lambda u,v: cmp(degree_within(G, u, within),
                                degree_within(G, v, within)))


def second_dist_to_set(G, v, W):
  """
  For vertices adjacent to W, return the lenght of the second shortest
  path to W. Otherwise return G.order(). WARNING: Not very efficient.
  """
  if v in W:
    return 0
  # Neighbors in W 
  nw = set(G.neighbors(v)).intersection(set(W))
  if len(nw) == 0:
    return G.order()
  if len(nw) >= 2:
    return 1
  G2 = Graph(G) # always undirected, copy
  G2.delete_edge(v, nw.pop())
  ds2 = G2.distance_all_pairs()
  dist2 = min([ds2[v][u] for u in W])
  assert dist2 >= 2
  return dist2


def order_max_adjacent(G, ordered = None, priorities = ['within', 'dist2', 'degree']):
  """
  Order the vertices of G trying to maximize the number of edges from each vertex
  to their predecessors. Starts with empty sequence and always appends a vertex maximizing
  the priorities given by `priorities` (W denotes the set of vertices already in the
  sequence or in `ordered`):

  * 'within' prefers vertices with most neighbors in W.
  * 'dist2' sorts first vertices adjacent to W by the second shortest path to W, then 
    all vertices nonadjacent to W. Might be slower to compute.
  * 'degree' prefers vertices with larger degree in `G`.
  * 'random' prefers random vertices (randomizes the order)

  If given, `ordered` is taken to be an already preordered sequence of vertices.
  Only the vertices outside `ordered` are then ordered and returned.

  Returns a list of vertices.
  """

  ordered = set(ordered or [])
  assert ordered.issubset(set(G.vertices()))

  if len(ordered) == G.order():
    return []

  vs = [v for v in G.vertices() if v not in ordered]
  for p in reversed(priorities):
    if p == 'within':
      sort_vertices_by_degree(G, vs, within=ordered)
    elif p == 'degree':
      sort_vertices_by_degree(G, vs)
    elif p == 'ranfom':
      shuffle(vs)
    elif p == 'dist2':
      # Cache the second distances
      d2 = dict([(v, second_dist_to_set(G, v, ordered)) for v in vs])
      vs.sort(cmp=lambda u,v: cmp(d2[u], d2[v]), reverse=True)
    else:
      raise Exception('Unknown priority "%s"', p)

  v = vs[-1]
  return [v] + order_max_adjacent(G, ordered=ordered.union([v]), priorities=priorities)


### Graph <-> CGraph conversion


def cgraph_to_graph(G, cls=Graph):
  """
  Convert a CGraph (without the fancy envelope) to a "normal"
  Sage Graph (or other class given by `cls`).
  """
  G2 = cls()
  G2.add_vertices(G.verts())
  for v in G.verts():
    for u in G.out_neighbors(v):
      G2.add_edge(v, u)
  return G2


def graph_to_cgraph(G, vmap=None, imap=None, graphtype=DenseGraph):
  """
  Convert given Graph to a CGraph (DenseGraph by default) on vertices [0 .. (order-1)]
  If provided, vmap is filled with {vertex: index}.
  If provided, imap is filled with {index: vertex}.
  Returns the CGraph
  """
  if vmap is None:
    vmap = {}
  if imap is None:
    imap = {}
  vs = list(G.vertices())

  CG = graphtype(G.order())
  for vi in range(len(vs)):
    v = vs[vi]
    vmap[v] = vi
    imap[vi] = v

  for e in G.to_directed().edges():
    CG.add_arc(vmap[e[0]], vmap[e[1]])
  return CG


### Homomorphism images and cores


def find_hom_image(G, candidates=None):
  """
  Tries to find a smaller hom-image of G. Repeatedly sets H:=G\{v} for v in candidates.
  Returns an induced subgraph on a strictly smaller hom-image, or None if G is a core.
  candidates defaults to V(G). If G is e.g. vertex-transitive, you may wish to use just
  one vertex.
  """
  if candidates is None:
    candidates = G.vertices()
  for v0 in candidates:
    H = G.copy()
    H.delete_vertex(v0)
    maps = extend_hom(G, H, limit=1)
    if len(maps) >= 1:
      VC = maps[0].values()
      C = G.subgraph(VC)
      return C
  return None


def check_subcore(G, candidates=None):
  """
  Tries to find a core, returning a smaller core or None if G is already core.
  Tries to remove every vertex of tryremove (one at a time), defauts to all.
  """
  if candidates is None:
    candidates = G.vertices()
  #logging.info(tryremove)
  for v0 in candidates:
    H = G.copy()
    H.delete_vertex(v0)
    maps = extend_hom(G, H, limit=1)
    if len(maps) >= 1:
      VC = maps[0].values()
      C = G.subgraph(VC)
      return C
  return None


def find_core(G, vertex_trans=False):
  """
  Finds a core of a Graph. If G is vertex-transitive, the process is much faster.
  Returns a core (a subgraph of G) or G itself.
  Complete graphs are treated specially.
  """
  if G.size() == (G.order() * (G.order()-1)) / 2:
    return G
  Gold = G
  rem0 = None
  if vertex_trans:
    rem0 = G.vertices()[:1]
  G = check_subcore(G, candidates=rem0)
  while G:
    Gold = G
    G = find_hom_image(G)
  return Gold


### Verifying homomorphisms


def is_hom(G, H, hom):
  """
  Check whether hom is a homomorphism from G to H.
  Works for both directed and undirected.
  """
  assert set(G.vertices()) == set(hom.keys())
  assert set(H.vertices()).issuperset(set(hom.values()))
  for e in G.edges():
    u, v = e[0], e[1]
    if not H.has_edge(hom[u], hom[v]):
      return False
  return True


### Module unit testing


def raises(f):
  "Internal helper for testing."
  try:
    f()
  except:
    return True
  return False


def test():
  "Run unit tests for the module."
  logging.getLogger().setLevel(logging.DEBUG)

  from sage.graphs.graph_generators import graphs
  K2 = graphs.CompleteGraph(2)
  K4 = graphs.CompleteGraph(4)
  K24 = graphs.CompleteGraph(24)
  C4 = graphs.CycleGraph(4)
  C16 = graphs.CycleGraph(16)

  # extend_hom
  assert len(extend_hom(K4, K4, limit=100)) == 24
  assert len(extend_hom(K2, K4, partmap={0:0}, limit=10)) == 3
  assert len(extend_hom(C16, K2, limit=10)) == 2
  assert len(extend_hom(C16, K2, partmap={0:0, 2:1}, limit=10)) == 0

  # extend_hom giving partial mappings
  assert len(extend_hom(C16, C4, partmap={0:0}, order=[0,2,1], limit=100)) == 4
  assert len(extend_hom(C16, C4, partmap={8:0}, order=[0], limit=100)) == 4
  assert extend_hom(C16, C4, partmap={4:1}, order=[], limit=100) == [{4:1}]
  assert extend_hom(C16, C4, order=[], limit=100) == [{}]

  # extend_hom only counting
  assert extend_hom(K2, K4, partmap={0:0}, limit=10, onlycount=True) == 3
  # Following would probably segfault if actually allocating 1G entries
  assert extend_hom(K4, K4, limit=2**30, onlycount=True) == 24

  # extend_hom with check_automorphisms
  assert len(extend_hom(K4, K4, limit=100, check_automorphisms=1)) == 6
  assert len(extend_hom(K4, K4, limit=100, check_automorphisms=2)) == 2
  assert len(extend_hom(K4, K4, limit=100, check_automorphisms=3)) == 1
  assert len(extend_hom(K4, K4, limit=100, check_automorphisms=4)) == 1
  assert len(extend_hom(K2, K4, partmap={0:0}, limit=10, check_automorphisms=1)) == 1
  assert len(extend_hom(C16, K2, limit=10, check_automorphisms=1)) == 1
  assert len(extend_hom(C16, K2, partmap={0:0, 2:1}, limit=10, check_automorphisms=2)) == 0
  # This might be slow
  assert extend_hom(K24, K24, limit=10, check_automorphisms=42, onlycount=True) == 1

  # find_hom_image
  K4c = find_hom_image(K4.disjoint_union(K4))
  assert K4.is_isomorphic(K4c)

  # find_core
  K4d = find_core(K4.disjoint_union(K4).disjoint_union(K4))
  assert K4.is_isomorphic(K4d)

  # is_hom
  for hom in extend_hom(C16, K2, limit=10):
    assert is_hom(C16, K2, hom)
  assert raises(lambda: is_hom(K2, K2, {0:0})) # not all vertices mapped
  assert raises(lambda: is_hom(K2, K2, {0:4})) # some vertices outside K2
  assert not is_hom(K2, K2, {0:0, 1:0}) # creates a loop

  logging.info("All tests passed.")


def benchmark():
  "Run few simple benchmarks for extend_hom()."
  logging.getLogger().setLevel(logging.DEBUG)

  from sage.graphs.graph_generators import graphs
  from sage.misc.sage_timeit import sage_timeit

  K6 = graphs.CompleteGraph(6)
  K2 = graphs.CompleteGraph(2)
  K7_6 = graphs.CompleteBipartiteGraph(7, 6)
  C16 = graphs.CycleGraph(16)
  C13 = graphs.CycleGraph(13)
  C5 = graphs.CycleGraph(5)
  eh = extend_hom

  logging.info("All homs K6 -> K6: %s", sage_timeit(
    "assert len(eh(K6, K6, limit=10000)) == 720", locals(), preparse=False))

  logging.info("All homs K7,6 -> K2: %s", sage_timeit(
    "assert len(eh(K7_6, K2, limit=10000)) == 2", locals(), preparse=False))

  logging.info("All homs C13 -> C5: %s", sage_timeit(
    "assert len(eh(C13, C5, limit=10000)) == 7150", locals(), preparse=False))

