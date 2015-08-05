import networkx as nx
import homsearch

G1 = nx.Graph([
    ('A','B'),
    ('A','C'),
    ('A','D'),
    ('A','E'),
    ('B','C'),
    ('C','D'),
    ('D','E'),
    ])

### Retracts

assert homsearch.find_retracts(G1, only_count=True) == 6

R1 = homsearch.find_retracts(G1, only_count=False)
assert R1[0]['A'] == 'A'

assert homsearch.find_retracts(G1, only_count=True, partmap={'A':'B'}) == 0

assert homsearch.find_retracts(G1, only_count=True, partmap={'B':'B'}) == 3

### Homomorphisms

assert homsearch.find_homomorphisms(G1, G1, only_count=True) == 36

R2 = homsearch.find_homomorphisms(G1, G1, only_count=False, partmap={'A':'E', 'E':'D', 'D':'A'})
assert (len(R2) == 1) and (R2[0]['C'] == 'D') and (R2[0]['B'] == 'A')

