"""
Cython routines for finding graph homomorphisms.
Do not use directly, only with the python module.
"""

from sage.graphs.base.dense_graph import DenseGraph


from sage.graphs.base.dense_graph cimport CGraph

cdef extern from "string.h":
  void *memcpy(void *dest, void *src, size_t n)

cdef extern from "alloca.h":
  void *alloca(size_t size)


def extend_hom(G, H, partmap = None, order = None, limit = 1,
               Ggraphtype=DenseGraph, Hgraphtype=DenseGraph,
               check_automorphisms=0, onlycount=False):
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
    check_automorphisms    If >0, check H for authomorphisms (after coloring
              the already targeted vertices) and only try one vertex per orbit.
              check_automorphisms is decreased with every recursion level and
              set to 0 on detection of all singleton orbits.
    onlycount  do not store all mappings, only return the number (still up to limit)

  Returns:
    list of mappings {v: f(v)} or their number (with `onlycout` set)
  """

  from homomorphisms import graph_to_cgraph, order_max_adjacent

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

  cdef int** resmaps_c = NULL
  if not onlycount:
    resmaps_c = <int **>alloca(sizeof(int*) * limit)
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

  r = extend_hom_c(CG, CH, partmap_c, order_c, len(order), resmaps_c, limit,
                   check_automorphisms)

  if onlycount:
    return r

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
                      int** resmaps, int limit,
                      int check_automorphisms):
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
    check_automorphisms    If >0, check H for authomorphisms (after coloring
              the already targeted vertices) and only try one vertex per orbit.
              check_automorphisms is decreased with every recursion level and
              set to 0 on detection of all singleton orbits.
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

  # array of vertices of H to try (in that order)
  cdef int *hverts = <int *>alloca(sizeof(int) * m)
  for i in range(m):
    hverts[i] = i
  cdef int hverts_c = m

  # H-automorphism elimination, only picks one vertex to try per
  # orbit of H. Potentially inefficient.
  if check_automorphisms > 0:
    from homomorphisms import cgraph_to_graph

    # construct vetrex partition isolating already used vertices of H
    Hused = set([partmap[Gv] for Gv in range(n)])
    Hparts = []
    Hfree = []
    for Hv in range(m):
      if Hv in Hused:
        Hparts.append([Hv])
      else:
        Hfree.append(Hv)
    Hparts.append(Hfree)

    # automorphism group and orbits of H
    HO = cgraph_to_graph(H)
    Horbits = HO.automorphism_group(partition=Hparts, orbits=True)[1]

    # set the vertices to try
    hverts_c = len(Horbits)
    for Hi in range(hverts_c):
      hverts[Hi] = Horbits[Hi][0]

    # condition to stop trying to find authomorphisms
    if hverts_c == m:
      check_automorphisms = 0

  # try all possibilities for fv=f(v)
  cdef int numres = 0
  cdef int fvi, fv, ok, u, fu, r
  for fvi in range(hverts_c):
    fv = hverts[fvi]
    ok = 1
    partmap2[v] = fv
    # check homomorphism corectness between v and every already mapped vertex
    # Treat as a dense graph for simplicity and speed 
    for u in range(n):
      fu = partmap2[u]
      if fu >= 0:
        if G.has_arc_unsafe(u, v) and (H.has_arc_unsafe(fu, fv) == 0):
          ok = 0
        if G.has_arc_unsafe(v, u) and (H.has_arc_unsafe(fv, fu) == 0):
          ok = 0
    # Recurse
    if ok:
      r = extend_hom_c(G, H, partmap2, tomap + 1, tomap_c - 1, resmaps, limit,
                       check_automorphisms - 1)
      numres += r
      limit -= r
      if resmaps != NULL:
        resmaps += r
  return numres

