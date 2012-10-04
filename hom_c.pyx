
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


def order_max_adjacent(G, ordered = None):
  if ordered is None:
    ordered = set([])
  ordered = set(ordered)
  assert ordered.issubset(set(G.vertices()))

  if len(ordered) == G.order():
    return []

  vmax = None
  adjmax = -1
  for v in G.vertices():
    if v not in ordered:
      adj = len([u for u in G.neighbors(v) if u in ordered])
      if adj > adjmax:
        adjmax = adj
        vmax = v
  assert adjmax >= 0

  return [vmax] + order_max_adjacent(G, ordered.union([vmax]))


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
