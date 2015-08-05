# Copyright (c) 2015 Tomas Gavenciak <gavento@ucw.cz>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to
# deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
# IN THE SOFTWARE.

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

