import collections
import logging
import random

from sage.parallel.decorate import parallel

import homomorphisms as hom


log = logging.getLogger("hom")


@parallel() # Guess no of CPUs
def parallel_extend_hom_partmaps(G, H, partmap, kwargs):
  log.debug("Running extend_hom with %s" % (partmap))
  res = hom.extend_hom(G, H, partmap=partmap, **kwargs)
  log.debug("Done    extend_hom with %s" % (partmap))
  return res


def parallel_extend_hom(G, H, partmap=None, branchdepth=3, onlycount=False,
                        order=None, limit=1, branchlimit=10000,
                        check_automorphisms=0, shuffle_jobs=False,
                        **kwargs):
  """
  Find homomorphism extensions in parallel using `extend_hom`.
  
  Starting with a given `partmap` (empty by default), extends the partial
  mappings to more and more vertices with parallel `extend_hom` jobs.
  The progression of extensions is given by `branchdepth` (either a list or
  a single number), additionally there is always an iteration extending the
  mappings to all vertices from `order` (all by default).
  Only up to `branchlimit` partial mappngs are generated in each iteration,
  so with high `branchdepth`, not all possibilities may be examined.
  The partial mappings are generated with respect to `order` and
  `check_automorphisms`, `check_automorphisms` is adjusted for each iteration
  to simulate `extend_hom` behavior.
  
  All other parameters are as for `extend_hom`.

  Note that each final `extend_hom` job is given the same `limit` (as there is
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
    order = hom.order_max_adjacent(G, partmap.keys())
  if not isinstance(branchdepth, collections.Iterable):
    branchdepth = [branchdepth]
  else:
    branchdepth = list(branchdepth)
  branchdepth.append(len(order))

  # Iteratively construct partmaps and then full map
  partmaps = [partmap]
  for bdi in range(len(branchdepth)):
    bd = branchdepth[bdi]
    finalrun = (bdi == len(branchdepth) - 1)

    # lenght of current partial maps (all the same) 
    curlen = len(partmaps[0])
    mapped = partmaps[0].keys()
    assert all([len(i) == curlen for i in partmaps])
    assert all([i.keys() == mapped for i in partmaps])

    # adjust branch depth to be at most order
    bd = min(bd, len(order))
    log.info("Finding partial maps of depth %d (of %d) from %d partial maps of depth %d." %
             (bd, len(order), len(partmaps), curlen))

    # Order restricted to unmapped vertices 
    curorder = [x for x in order if x not in mapped]

    # Maximum number of maps to look for (upper bound)
    maxmaps = max(len(H.vertices())**bd, len(partmaps)) + 1

    # Construct all job arguments
    args = []
    for pm in partmaps:
      kwcopy = kwargs.copy()
      kwcopy.update(onlycount=False, order=order[:bd], limit=min(maxmaps, branchlimit))
      # last iteration?
      if finalrun:
        kwcopy.update(onlycount=onlycount, limit=limit)
      # proper depth of automorphism checking
      kwcopy.update(check_automorphisms = max(check_automorphisms + len(partmap) - curlen, 0))
      args.append((G, H, pm, kwcopy))

    # Experimental - shuffle the job order
    if shuffle_jobs:
      random.shuffle(args)

    # Run the jobs
    rescount = 0
    res = []
    done = 0
    jobs = parallel_extend_hom_partmaps(args)

    log.info(" %4d of %d done, %d results" % (done, len(args), rescount))

    # Monitor the jobs as they are finished
    for k,v in jobs:
      if v == 'NO DATA':
        raise Exception("Error in worker, no data returned.")
      if finalrun and onlycount:
        rescount += v
      else:
        rescount += len(v)
        res.extend(v)
      done += 1

      nth = max(len(args) / (100 if finalrun else 20), 1)
      allupto = 5
      if done <= allupto or (done % nth == 0) or rescount >= limit:
        log.info(" %4d of %d done, %d results" % (done, len(args), rescount))
      if done == allupto and nth > 1:
        log.info(" [from now on reporting only every %d-th]" % (nth))

      # Watch for end conditions on final run
      if finalrun:
        if rescount >= limit or done >= len(args):
          del jobs
          if onlycount:
            return min(rescount, limit)
          else:
            return res[:limit]

    assert not finalrun
    # For partial maps (final run should never get here)
    assert rescount == len(res)

    if rescount == 0:
      if onlycount:
        return 0
      else:
        return []

    assert rescount < maxmaps
    if rescount >= branchlimit:
      log.error("Branchlimit reached - not all homomorphisms will be searched. Consider decreasing branchdepth.")
    partmaps = res

  # This should never be reached
  assert False

def test():
  "Some unit tests"

  from sage.graphs.graph_generators import graphs
  K2 = graphs.CompleteGraph(2)
  K4 = graphs.CompleteGraph(4)
  K24 = graphs.CompleteGraph(24)
  C4 = graphs.CycleGraph(4)
  C16 = graphs.CycleGraph(16)

  # parallel_extend_hom
  assert len(parallel_extend_hom(K4, K4, limit=100)) == 24
  assert parallel_extend_hom(K2, K4, partmap={0:0}, limit=10, onlycount=True) == 3
  assert len(parallel_extend_hom(C16, K2, limit=10)) == 2
  assert len(parallel_extend_hom(C16, K2, partmap={0:0, 2:1}, limit=10)) == 0

  # parallel_extend_hom giving partial mappings
  assert len(parallel_extend_hom(C16, C4, partmap={0:0}, order=[0,2,1], limit=100)) == 4
  assert len(parallel_extend_hom(C16, C4, partmap={8:0}, order=[0], limit=100)) == 4
  assert parallel_extend_hom(C16, C4, partmap={4:1}, order=[], limit=100) == [{4:1}]
  assert parallel_extend_hom(C16, C4, order=[], limit=100) == [{}]

  # differen branching parameters
  assert parallel_extend_hom(C4, K24, branchdepth=[0,0,1], limit=1000000, onlycount=True, check_automorphisms=1) == 11661
  assert parallel_extend_hom(C4, K24, branchdepth=[0,1,2,3], limit=1000000, onlycount=True, check_automorphisms=2) == 507
  assert (sorted(parallel_extend_hom(C4, K24, branchdepth=[1,0,2,3,2], limit=10000, check_automorphisms=14)) ==
      [{0: 1, 1: 0, 2: 1, 3: 0}, {0: 1, 1: 2, 2: 1, 3: 0}, {0: 2, 1: 0, 2: 1, 3: 0}, {0: 3, 1: 2, 2: 1, 3: 0}])

  log.info("All tests passed.")

