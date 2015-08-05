def ricardo_graph(n):
  r""" 
    Returns a graph I learned about from Ricardo Strausz, 
    but they are already mentioned in a paper 
    Aracho, Bracho, Neumann von Lara: ??? Tight hypergraphs 
  """ 

  p = 2*n+1
  if not p.is_prime():
#    print "2*n+1 = %d should be a prime"%(2*n+1) 
    return None 


  Z = IntegerModRing(p) 
  G = Graph() 

  for x in [1..n]:
    for y in [x+1..n]:
      x = Z(x) 
      y = Z(y) 
      for z in third(x,y):
        G.add_vertex(frozenset((x,y,z))) 

  for v in G.vertices():
    for x in v:
      y,z = v.difference({x}) 
      for w in set(third(y,z)): 
        if w != x:
          neigh = frozenset((w,y,z))
          G.add_edge(v,neigh) 

  return G 


def myabs(x):
  r"""
    x is element of some Z_p
    returns x or -x, whichever is smaller.       
  """ 

  return min(x,-x) 

def third(x,y):
  ## does not work for x=0 
  x = myabs(x) 
  y = myabs(y) 
  Z = x.parent() 
  if x==y:
    return [] 
  if x==0 or y==0:
    x = x+y
    return [ myabs(2*x), myabs(x/2) ] 
  L = [] 
  for z in [x-y, x+y]: 
    z = myabs(z) 
    if z==x or z==y: 
      L.append(0)
    else:
      L.append(z) 
  return L 


