
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


def extend_hom(G, H, partmap = None, reslimit = 1):
  """
  Recursive routine to extend partial homomorphism G->H given by partmap to a
  full homomorphism. Finds up to reslimit solutions. Works correctly on digraphs.
  If a partial map is provided, homomorphism correctness is not checked on set vertices.

  Parameters:
    G         source Graph
    H         target Graph
              Note that the vertices should be hashable values.
    partmap   partial map {v: f(v)}, empty by default
    reslimit  maximum number of homomorphisms to look for, default 1

  Returns:
    list of mappings {v: f(v)}
  """

  if partmap is None:
    partmap = {}

  # Preprocess to CGraphs
  Gimap = {}
  Gvmap = {}
  CG = graph_to_cgraph(G, vmap=Gvmap, imap=Gimap)
  Himap = {}
  Hvmap = {}
  CH = graph_to_cgraph(H, vmap=Hvmap, imap=Himap)

  # General init
  n = len(CG.verts())

  cdef int** resmaps_c = <int **>alloca(sizeof(int*) * reslimit)
  for i in range(reslimit):
    resmaps_c[i] = <int *>alloca(sizeof(int) * n)

  cdef int *partmap_c = <int *>alloca(sizeof(int) * n)
  for v in G.vertices():
    if v in partmap:
      partmap_c[Gvmap[v]] = Gvmap[partmap[v]]
    else:
      partmap_c[Gvmap[v]] = -1

  r = extend_hom_c(CG, CH, partmap_c, NULL, -1, NULL, resmaps_c, reslimit)
  res = []
  for i in range(r):
    new = {}
    for j in range(n):
      # Convert vertex mapping back to back
      new[Gimap[j]] = Himap[resmaps_c[i][j]]
    res.append(new)

  return res

cdef int extend_hom_c(CGraph G, CGraph H, int partmap[],
                      int tomap[], int tomap_c, int adjdeg[],
                      int** resmaps, int reslimit):
  """
  Recursive routine to extend partial homomorphism G->H given by partmap to a
  full homomorphism. Finds up to reslimit solutions, storing them in resmaps.
  Works correctly on digraphs.

  Parameters:
    G         source CGraph
    H         target CGraph
    partmap   partial map, partmap[v]=f(v), must be -1 for unmapped vertices
              if NULL, created as empty map
    tomap     list of unmapped vertices (in positions 0..tomap_c-1)
              if NULL, created according to partmap
    tomap_c   number of vertices to map, ignored if tomap==NULL
    resmaps   array of int[n] pointers to store complete maps,
              size must be at least reslimit
              if NULL, no results are written (just counted)
    reslimit  maximum number of homomorphisms to look for

  Return value: number of homomorphisms found and returned in resmaps
  """

  cdef int n = G.num_verts
  cdef int m = H.num_verts
  cdef int i, j

  if reslimit <= 0:
    return 0

  if tomap == NULL:
    tomap = <int *>alloca(sizeof(int) * n)
    tomap_c = 0
    for i in range(n):
      if partmap[i] < 0:
        tomap[tomap_c] = i
        tomap_c += 1

  if partmap == NULL:
    partmap = <int *>alloca(sizeof(int) * n)
    for i in range(n):
      partmap[i] = -1

  if adjdeg == NULL:
    adjdeg = <int *>alloca(sizeof(int) * n)
    for i in range(n):
      adjdeg[i] = 0
    for i in range(n):
      if partmap[i] >= 0:
        for j in range(n):
          if G.has_arc(i, j):
            adjdeg[j] += 1
          if G.has_arc(j, i):
            adjdeg[j] += 1

  if tomap_c == 0:
    if resmaps != NULL:
      memcpy(resmaps[0], partmap, n * sizeof(int))
    return 1


  # pick a node index to map
  # TODO: be smarter ;-)
  vi = 0
  cdef int deg, vi_deg = 0
  for i in range(tomap_c):
    deg = 0
    for j in range(n):
      if partmap[j] >= 0:
        if G.has_arc(tomap[i], j):
          deg += 1
        if G.has_arc(j, tomap[i]):
          deg += 1
    #print ("%d,%d" % (deg, adjdeg[tomap[i]])),
    assert deg == adjdeg[tomap[i]]
    if deg > vi_deg:
      vi_deg = deg
      vi = i

  v = tomap[vi]
  print v
  # copy tomap to tomap2, remove v from tomap2
  cdef int *tomap2 = <int *>alloca(sizeof(int) * n)
  memcpy(tomap2, tomap, sizeof(int) * n)
  tomap2[vi] = tomap2[tomap_c-1]
  # copy partmap to partmap2
  cdef int *partmap2 = <int *>alloca(sizeof(int) * n)
  memcpy(partmap2, partmap, sizeof(int) * n)
  # copy of adjdeg
  cdef int *adjdeg2 = <int *>alloca(sizeof(int) * n)
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
        if G.has_arc(u, v) and (H.has_arc(fu, fv) == 0):
          ok = 0
        if G.has_arc(v, u) and (H.has_arc(fv, fu) == 0):
          ok = 0
    # Recurse
    if ok:
      # update candidate degrees
      memcpy(adjdeg2, adjdeg, sizeof(int) * n)
      for j in range(n):
        if partmap2[j] < 0:
          if G.has_arc(j, v):
            adjdeg2[j] += 1
          if G.has_arc(v, j):
            adjdeg2[j] += 1

      r = extend_hom_c(G, H, partmap2, tomap2, tomap_c-1, adjdeg2, resmaps, reslimit)
      numres += r
      reslimit -= r
      if resmaps != NULL:
        resmaps += r
  return numres
