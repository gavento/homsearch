import homomorphisms as hom
import logging as log
from sage.parallel.decorate import parallel

@parallel() # Guess no of CPUs
def parallel_extend_hom_partmaps(G, H, partmap, kwargs):
  log.info("Running extend_hom with %s"%(partmap))
  res = hom.extend_hom(G, H, partmap=partmap, **kwargs)
  log.info("Done    extend_hom with %s"%(partmap))
  return res


def parallel_extend_hom(G, H, partmap=None, branchdepth=3, onlycount=False, order=None, limit=1, branchlimit=100000, **kwargs):
  if partmap is None:
    partmap = {}
  if order is None:
    # Only nonmapped vertices
    order = [x for x in hom.order_max_adjacent(G, partmap.keys()) if x not in partmap]
  branchdepth = min(branchdepth, len(order))

  partmaps = hom.extend_hom(G, H, partmap=partmap, order=order[:branchdepth],
                            onlycount=False, limit=min(len(H.vertices())**branchdepth, branchlimit), **kwargs)
  if len(partmaps) >= branchlimit:
    log.error("Branchlimit reached - not all homomorphisms will be searched. Consider decreasing branchdepth.")

  print("Found %d partmaps" % len(partmaps))
  args = []
  for pm in partmaps:
    kwcopy = kwargs.copy()
    kwcopy.update(onlycount=onlycount, order=order, limit=limit)
    if 'check_automorphisms' in kwcopy:
      kwcopy['check_automorphisms'] = max(kwcopy['check_automorphisms'] - branchdepth, 0)
    args.append((G, H, pm, kwcopy))

  res = list(parallel_extend_hom_partmaps(args))
  res2 = [b for a,b in res]

  if onlycount:
    return min(sum(res2), limit)
  else:
    return flatten(res2)[:limit]

