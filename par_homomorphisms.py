import homomorphisms as hom
import logging as log
from sage.parallel.decorate import parallel

@parallel() # Guess no of CPUs
def parallel_extend_hom_partmaps(G, H, partmap, kwargs):
  log.debug("Running extend_hom with %s"%(partmap))
  res = hom.extend_hom(G, H, partmap=partmap, **kwargs)
  log.debug("Done    extend_hom with %s"%(partmap))
  return res


def parallel_extend_hom(G, H, partmap=None, branchdepth=3, onlycount=False, order=None, limit=1, branchlimit=10000, **kwargs):
  """
  Find homomorphism extensions in parallel using `extend_hom`.
  
  First computes all partial maps of first `branchdepth` vertices and then
  extends each mapping using `extend_hom` as a parallel task. Only up to
  `branchlimit` partial mappngs are generated, so with high `branchdepth`, not
  all may be examined. The partial mappings are generated with respect to
  `order` and `check_automorphisms`.
  
  All other parameters are as for `extend_hom`.

  Note that each branched `extend_hom` is given the same `limit` (as there is
  no guessing which branch will have enough results) and therefore
  the number of results computed may be substantially higher (but not more
  than `limit` is returned). As soon as the goal `limit` is met, the remaining
  jobs are terminated.
  
  Uses Sage `@parallel` with CPU numeber autodetection and forking.
  """

  if partmap is None:
    partmap = {}
  if order is None:
    # Only nonmapped vertices
    order = [x for x in hom.order_max_adjacent(G, partmap.keys()) if x not in partmap]
  branchdepth = min(branchdepth, len(order))

  maxmaps = len(H.vertices())**branchdepth + 1
  partmaps = hom.extend_hom(G, H, partmap=partmap, order=order[:branchdepth],
                            onlycount=False, limit=min(maxmaps, branchlimit), **kwargs)
  assert len(partmaps) < maxmaps
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

  rescount = 0
  res = []
  done = 0

  jobs = parallel_extend_hom_partmaps(args)

  log.info("%4d of %d done, %d results", done, len(args), rescount)

  for k,v in jobs:
    if onlycount:
      rescount += v
    else:
      rescount += len(v)
      res.extend(v)
    done += 1

    nth = min(len(args)/100, 1)
    if done <= 42 or (done % nth == 0) or rescount >= limit:
      log.info("%4d of %d done, %d results", done, len(args), rescount)
    if done == 42:
      log.info("[from now on reporting only every %d-th]", nth)

    if rescount >= limit:
      del jobs
      break

  if onlycount:
    return min(rescount, limit)
  else:
    return res[:limit]

