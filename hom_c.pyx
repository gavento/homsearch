from random import shuffle
from sage.graphs.graph import Graph
from sage.graphs.base.dense_graph import DenseGraph
from sage.graphs.base.dense_graph cimport CGraph

cdef extern from "string.h":
  void *memcpy(void *dest, void *src, size_t n)

cdef extern from "alloca.h":
  void *alloca(size_t size)


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


def order_max_adjacent(G, ordered = None, priorities=['within', 'dist2', 'degree']):
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
  return [v] + order_max_adjacent(G, ordered.union([v]))


def extend_hom(G, H, partmap = None, order = None, limit = 1,
               Ggraphtype=DenseGraph, Hgraphtype=DenseGraph):
  """
  Recursive routine to extend partial homomorphism G->H given by partmap to a
  full homomorphism. Finds up to limit solutions. Works correctly on digraphs.
  If a partial map is provided, homomorphism correctness is not checked on set vertices.

  Parameters:
    G         Source Graph.
    H         Target Graph.
              Note that the vertices of both G and H should be hashable values.
    partmap   Partial map {v: f(v)}, empty by default.
    order     Sequence of all the unmapped vertices of G giving the order for
              the backtracking. Already mapped vertices are ignored.
              Defaults to heuristic order_max_adjacent().
    limit     Maximum number of homomorphisms to look for and return, default 1.
    Ggraphtype, Hgraphtype     The CGraph class to use for C representation of G
              and H. DenseGraph by default.

  Returns:
    list of mappings {v: f(v)}
  """

  if partmap is None:
    partmap = {}

  if order is None:
    order = order_max_adjacent(G, partmap.keys())
  order = [v for v in order if v not in partmap.keys()]
  assert set(order).union(set(partmap.keys())) == set(G.vertices())

  # Preprocess to CGraphs
  Gimap = {}
  Gvmap = {}
  CG = graph_to_cgraph(G, vmap=Gvmap, imap=Gimap, graphtype=Ggraphtype)
  Himap = {}
  Hvmap = {}
  CH = graph_to_cgraph(H, vmap=Hvmap, imap=Himap, graphtype=Hgraphtype)

  # General init
  n = len(CG.verts())

  cdef int** resmaps_c = <int **>alloca(sizeof(int*) * limit)
  for i in range(limit):
    resmaps_c[i] = <int *>alloca(sizeof(int) * n)

  cdef int *partmap_c = <int *>alloca(sizeof(int) * n)
  for v in G.vertices():
    if v in partmap:
      partmap_c[Gvmap[v]] = Hvmap[partmap[v]]
    else:
      partmap_c[Gvmap[v]] = -1

  cdef int *order_c = <int *>alloca(sizeof(int) * n)
  for i in range(len(order)):
    order_c[i] = Gvmap[order[i]]

  r = extend_hom_c(CG, CH, partmap_c, order_c, len(order), resmaps_c, limit)
  res = []
  for i in range(r):
    new = {}
    for j in range(n):
      # Convert vertex mapping back to G and H
      new[Gimap[j]] = Himap[resmaps_c[i][j]]
    res.append(new)

  return res

cdef int extend_hom_c(CGraph G, CGraph H, int partmap[],
                      int tomap[], int tomap_c,
                      int** resmaps, int limit):
  """
  Recursive routine to extend partial homomorphism G->H given by partmap to a
  full homomorphism. Finds up to limit solutions, storing them in resmaps.
  Works correctly on digraphs.

  Parameters:
    G         Source CGraph.
    H         Target CGraph.
    partmap   Partial map, partmap[v]=f(v), must be -1 for unmapped vertices.
              If NULL, created as empty map.
    tomap     List of unmapped vertices in the order to be assigned by the
              backtracking routine. Must not contain vertices already in partmap.
              If NULL, use the numeric order on the unmapped vertices of G.
    tomap_c   Number of vertices to map, ignored if tomap==NULL.
    resmaps   Array of int[n] pointers to store complete maps,
              size must be at least limit.
              If NULL, no results are written (just counted).
    limit     Maximum number of homomorphisms to look for.

  Return value: number of homomorphisms found and returned in resmaps
  """

  cdef int n = G.num_verts
  cdef int m = H.num_verts
  cdef int i, j

  if limit <= 0:
    return 0

  if partmap == NULL:
    partmap = <int *>alloca(sizeof(int) * n)
    for i in range(n):
      partmap[i] = -1

  if tomap == NULL:
    tomap = <int *>alloca(sizeof(int) * n)
    tomap_c = 0
    for i in range(n):
      if partmap[i] < 0:
        tomap[tomap_c] = i
        tomap_c += 1

  if tomap_c == 0:
    if resmaps != NULL:
      memcpy(resmaps[0], partmap, n * sizeof(int))
    return 1


  v = tomap[0]

  # copy partmap to partmap2
  cdef int *partmap2 = <int *>alloca(sizeof(int) * n)
  memcpy(partmap2, partmap, sizeof(int) * n)

  # try all possibilities for w=f(v)
  cdef int numres = 0
  cdef int fv, ok, u, fu, r
  for fv in range(m):
    ok = 1
    partmap2[v] = fv
    # check homomorphism corectness between v and every already mapped vertex
    # Treat as dense graph for simplicity and speed 
    for u in range(n):
      fu = partmap2[u]
      if fu >= 0:
        if G.has_arc_unsafe(u, v) and (H.has_arc_unsafe(fu, fv) == 0):
          ok = 0
        if G.has_arc_unsafe(v, u) and (H.has_arc_unsafe(fv, fu) == 0):
          ok = 0
    # Recurse
    if ok:
      r = extend_hom_c(G, H, partmap2, tomap + 1, tomap_c - 1, resmaps, limit)
      numres += r
      limit -= r
      if resmaps != NULL:
        resmaps += r
  return numres
