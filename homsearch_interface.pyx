# distutils: language = c++
# distutils: sources = homsearch_lib.cpp

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


"""
Cython routines for finding graph homomorphisms.
Do not use directly, only with the python module.
"""

from libcpp.vector cimport vector
from libcpp cimport bool
from cython.operator cimport dereference as deref

cdef extern from "homsearch_lib.h":

    cdef cppclass homsearch:
        # Graphs   
        vector[vector[int]] G, H

        # Result options
        long long int res_count
        long long int res_limit
        vector[vector[int]] res_list
        bool res_store

        # Other options
        bool retract_mode
        int max_depth

        # Search interface 
        void search_vector(vector[int] &f, int depth) nogil
        void search(int depth) nogil

    # helper to create right sized homsearch
    homsearch *new_homsearch(vector[vector[int]] &G, vector[vector[int]] &H,
            long long int res_limit, bool res_store, bool retract_mode_, int max_depth)


cdef class HomsearchInterface:
    "Mid-level interface to homsearch_lib working with graphs as neighbor lists on [0 .. n-1]"

    cdef homsearch *srch

    def __init__(self, G_adj, H_adj, res_limit, res_store, retract_mode, max_depth=-1):
        self.srch = new_homsearch(G_adj, H_adj, res_limit, res_store, retract_mode, max_depth)

    def search(self):
        "Search from an empty mapping"
        cdef homsearch *s = self.srch
        with nogil:
            s.search(0)

    def search_from(self, f):
        "Search from a given partial mapping"
        assert isinstance(f, list)
        assert int(len(f)) == self.srch.G.size()
        cdef vector[int] vf = f

        cdef homsearch *s = self.srch
        with nogil:
            s.search_vector(vf, 0)

    def result_list(self):
        "Return list of found maps (empty when res_store==False)"
        return self.srch.res_list

    def result_count(self):
        "Return number of found maps (up to res_limit)"
        return self.srch.res_count

    def __del__(self):
        del self.srch



