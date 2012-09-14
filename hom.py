import logging as log
import itertools

from sage.graphs.graph import Graph
from sage.misc.preparser import load
from sage.misc.misc import powerset

try:
  load("hom_c.pyx", locals())
except Exception as e:
  log.error("hom_c.pyx failed to load: %s", e)

def gen_cayley_z2(n, gens):
  G = Graph()
  ranges = [(0,1)] * n
  for v in itertools.product(*ranges):
    v = tuple(v)
    G.add_vertex(v)
    for g in gens:
      w = tuple([(x + y) % 2 for (x, y) in zip(v, g)])
      G.add_edge(v, w)
  return G

g16 = gen_cayley_z2(4, [
  (1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1),
  (1,1,0,0), (1,1,0,0), (1,0,1,0), (1,0,0,1), (0,1,1,0), (0,1,0,1), (0,0,1,1)])

def check_subcore(G, tryremove=None):
  """
  Tries to find a core, returning a smaller core or None if G is already core.
  Tries to remove every vertex of tryremove (one at a time), defauts to all.
  """
  if tryremove is None:
    tryremove = G.vertices()
  #log.info(tryremove)
  for v0 in tryremove:
    #log.debug("Removing %s", v0)
    H = G.copy()
    H.delete_vertex(v0)
    maps = extend_hom(G, H, reslimit=2)
    if len(maps) >= 1:
      VC = maps[0].values()
      C = G.subgraph(VC)
      return C
  return None

def find_core(G, vertex_trans=False):
  "Finds a core of a (optionally vertex-transitive) graph."
  Gold = G
  rem0 = None
  if vertex_trans:
    rem0 = G.vertices()[:1]
  G = check_subcore(G, tryremove=rem0)
  while G:
    Gold = G
    G = check_subcore(G)
  return Gold

def generate_cubes(n):
  ranges = [(0,1)] * n
  vs = list(map(tuple, itertools.product(*ranges)))
  basegens = [v for v in vs if sum(v) == 1]
  optgens = [v for v in vs if sum(v) >= 2]
  seen_graphs = set()
  gensets = list(powerset(optgens))
  c = 0
  for gset in gensets:
    c += 1
    if len(gset) >= len(gensets):
      continue
    gens = basegens + gset
    G = gen_cayley_z2(n, gens)

    canon = tuple(Graph(G).canonical_label().edges())
    if canon in seen_graphs:
      continue
    log.debug("%d of %d: Generators: %s" % (c, len(gensets), gens))
    seen_graphs.add(canon)

    yield (gset, G)


def try_cubes(n):
  res = []
  for gens, G in generate_cubes(n):
    C = find_core(G)
    msg = "Core: order %d, size %d, complete %s" %(C.order(), C.size(), C.size() == C.order() * (C.order() - 1) / 2)
    log.info(msg)
    res.append((gens, G, C, msg))
  return res
