# distutils: language = c++
# distutils: sources = homsearch-lib.cpp
"""
Cython routines for finding graph homomorphisms.
Do not use directly, only with the python module.
"""

from libcpp.vector cimport vector
from libcpp cimport bool
from cython.operator cimport dereference as deref

cdef extern from "homsearch-lib.h":

    cdef cppclass search_state32:
        vector[int] f
        search_state32(homsearch32 *search_, vector[int] *f_)

    cdef cppclass homsearch32:
        long long int res_count
        long long int res_limit
        vector[vector[int]] res_list
        bool res_store
        bool retract_mode
        homsearch32(vector[vector[int]] &G_, vector[vector[int]] &H_, long long int res_limit_, bool res_store_, bool retract_mode_, int max_depth_)
        void search(search_state32 &s, int depth)

cdef class Homsearch32:
    cdef homsearch32 *search32

    def __init__(self, G_adj, H_adj, res_limit, res_store, retract_mode, max_depth=-1):
        max_size = max(len(G_adj), len(H_adj))
        assert (max_size <= 32)
        self.search32 = new homsearch32(G_adj, H_adj, res_limit, res_store, retract_mode, max_depth)
        pass

    def search(self):
        cdef search_state32 *s = new search_state32(self.search32, NULL)
        self.search32.search(deref(s), 0)
        del s
        return self.search32.res_list

    def __del__(self):
        del self.search32

cdef class Homsearch512:
    pass

def Homsearch_create(G_adj, H_adj, res_limit, res_store, retract_mode, max_depth=-1):
    max_size = max(len(G_adj), len(H_adj))
    if max_size <= 32:
        return Homsearch32(G_adj, H_adj, res_limit, res_store, retract_mode, max_depth)
    if max_size <= 512:
        return Homsearch512(G_adj, H_adj, res_limit, res_store, retract_mode, max_depth)
    raise Exception("Graph too big for homcore lib (max size 512)")


