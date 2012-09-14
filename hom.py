#import sage
from sage.graphs.graph import Graph as SageGraph
from networkx import Graph
import logging as log
import itertools

log.getLogger().setLevel(log.DEBUG)

def extend_hom(G, H, partmap=None, VG=None, limit=1):
  if limit <= 0:
    return []
  if partmap is None:
    partmap = {}
  if VG is None:
    VG = G.nodes()
  if VG == []:
    return [partmap][:limit]
  res = []

  #if len(partmap) <= 2:
  #  log.debug(len(partmap))

  # Check preservation of adjacency of mapped neighbor edges
  def adjacency_OK(v, m):
    for w in G.neighbors(v):
      if w in m and not H.has_edge(m[v], m[w]):
        return False
    return True

  # Count number of mapped neighbors
  mapdegs=dict([(v, 0) for v in G.nodes()])
  for v in partmap:
    for n in G.neighbors(v):
      mapdegs[n] += 1
  # Find vertex with maximum num. of mapper neighbors
  maxd = 0
  maxv = VG[0]
  for v in VG:
    if mapdegs[v] > maxd:
      maxv = v
      maxd = mapdegs[v]

  v = maxv
  VG2 = [i for i in VG if i != v]
  partmap2 = dict(partmap)
  for u in H.nodes():
    partmap2[v] = u
    if adjacency_OK(v, partmap2):
      r = extend_hom(G, H, partmap2, VG2, limit)
      limit -= len(r)
      res.extend(r)
      if limit == 0:
        return res
  return res

def gen_cayley_z2(n, gens):
  G = Graph()
  ranges = [(0,1)] * n
  for v in itertools.product(*ranges):
    v = tuple(v)
    G.add_node(v)
    for g in gens:
      w = tuple([(x + y) % 2 for (x, y) in zip(v, g)])
      G.add_edge(v, w)
  return G

#g=Graph(z2_4.cayley_graph(generators=[a,b,c,d,a*b,a*c,a*d,b*c,b*d,c*d]))
g16 = gen_cayley_z2(4, [
  (1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1),
  (1,1,0,0), (1,1,0,0), (1,0,1,0), (1,0,0,1), (0,1,1,0), (0,1,0,1), (0,0,1,1)])

def check_subcore(G, tryremove=None):
  if tryremove is None:
    tryremove = G.nodes()
  #log.info(tryremove)
  for v0 in tryremove:
    #log.debug("Removing %s", v0)
    H = G.copy()
    H.remove_node(v0)
    maps = extend_hom(G, H, limit=2)
    if len(maps) >= 1:
      VC = maps[0].values()
      C = G.subgraph(VC)
      return C
  return None

def find_core(G):
  Gold = G
  G = check_subcore(G, tryremove=G.nodes()[:1]) ## Not necessaruly correc assumption ...
  while G:
    Gold = G
    G = check_subcore(G)
  return Gold

def try_generators(n):
  pass
  ranges = [(0,1)] * n
  vs = list(map(tuple, itertools.product(*ranges)))
  basegens = [v for v in vs if sum(v) == 1]
  optgens = [v for v in vs if sum(v) >= 2]
  from sage.misc.misc import powerset
  seen_graphs = set()
  res = []
  gensets = list(powerset(optgens))
  c = 0
  for gset in gensets:
    c += 1
    if len(gset) >= len(sensets):
      continue
    gens = basegens + gset
    G = gen_cayley_z2(n, gens)

    canon = tuple(SageGraph(G).canonical_label().edges())
    if canon in seen_graphs:
      continue
    log.info("%d of %d: Generators: %s" % (c, len(gensets), gens))
    seen_graphs.add(canon)

    C = find_core(G)
    msg = "Core: order %d, size %d, complete %s" %(C.order(), C.size(), C.size() == C.order() * (C.order() - 1) / 2)
    log.info(msg)
    res.append((gens, G, C, msg))
  return res
