"""
Experiments for hypothesis about cube-like graph cores.
"""

import logging as log
import hom

log.getLogger().setLevel(log.DEBUG)


def cube_cores(n):
  "Find the cores for all n-dim cubes. Returns [(gens, G, Core, message), ]"
  res = []
  for gens, G in generate_cubes(n):
    C = find_core(G)
    msg = "Core: order %d, size %d, complete %s" %(C.order(), C.size(), C.size() == C.order() * (C.order() - 1) / 2)
    log.info(msg)
    res.append((gens, G, C, msg))
  return res



def check_some_squash_or_core(G):
  """
  Check whether given n-dim cube-like graph is either core itself,
  or can be squashed in one direction. Returns True or False.

  Fails for G gen by TODO.
  """
  if len(G.vertices()) <= 1:
    return True
  for c in G.vertices():
    if sum(map(int,c)) > 0:
      h = hom.squash_cube(G, c)
      if h:
        log.debug("Squashable by %s", c)
        return True
  # No squash worked
  H = hom.find_core(G, vertex_trans=True)
  if H.order() == G.order():
    # G itself core
    return True
  return False


def check_prop_for_all(prop, n):
  """
  Check prop(G)->True/False for all n-dim cube-like graphs.
  For example: check_prop_for_all(check_some_squash_or_core, 4)
  Returns list of failed graphs as pairs (gens, G).
  """
  fails = []
  for gen, G in hom.nonisomorphic_cubes_Z2(n):
    if not prop(G):
      log.debug("Prop failed: G with gen %s", gen)
      fails.append((gen, G))
  return fails


def check_cores_are_cubelike(n):
  smaller_cores_canon = set()
  for i in range(1,n):
    for gen, G in hom.nonisomorphic_cubes_Z2(i):
      Gcanon = tuple(G.canonical_label().edges())
      smaller_cores_canon.add(Gcanon)

  for gen, G in hom.nonisomorphic_cubes_Z2(n):
    H = hom.find_core(G)
    if H.order() == G.order():
      print "G itself core, gen by", gen
      continue
    Hcanon = tuple(hom.Graph(H).canonical_label().edges())
    if Hcanon in smaller_cores_canon:
      continue
    print "Not a cube-like core", H, "of", G, "gen by", gen
    return False
  return True


def check_some_core_hom_is_by_subspace(vs, G, H=None, require_isomorphic=True, limit=10000):
  """
  Check if there is a retract from G to H that "contracts" a subspace of vs.
  See hom.hom_is_by_subspace. If not given, H is a core of G.
  Returns True/False.
  """
  if not H:
    H = hom.find_core(G)
  return any((hom.hom_is_by_subspace(G, vs, h) for h in hom.extend_hom(G, H, limit=limit)))


