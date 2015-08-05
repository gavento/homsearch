/*
 * Copyright (c) 2015 Tomas Gavenciak <gavento@ucw.cz>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#include "homsearch_lib.h"

///////////////////////////////////////////////////////////
// Helper to create the right instance of homsearch_impl<>

homsearch *new_homsearch(const vector<vector<int> > &G, const vector<vector<int> > &H,
              long long int res_limit, bool res_store, bool retract_mode, int max_depth)
{
    int max_size = max(G.size(), H.size());

    if (max_size <= 16)
        return new homsearch_impl<16>(G, H, res_limit, res_store, retract_mode, max_depth);

    if (max_size <= 32)
        return new homsearch_impl<32>(G, H, res_limit, res_store, retract_mode, max_depth);

    if (max_size <= 64)
        return new homsearch_impl<64>(G, H, res_limit, res_store, retract_mode, max_depth);

    if (max_size <= 128)
        return new homsearch_impl<128>(G, H, res_limit, res_store, retract_mode, max_depth);

    if (max_size <= 256)
        return new homsearch_impl<256>(G, H, res_limit, res_store, retract_mode, max_depth);

    if (max_size <= 1024)
        return new homsearch_impl<1024>(G, H, res_limit, res_store, retract_mode, max_depth);

    if (max_size <= 4096)
        return new homsearch_impl<4096>(G, H, res_limit, res_store, retract_mode, max_depth);

    throw logic_error("homsearch_impl not implemented for graphs larger than 4096");
}


