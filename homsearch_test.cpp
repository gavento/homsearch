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

#include <iostream>

#define push_array(G, A) { const vector<int> pav(A, A + (sizeof(A) / sizeof(int))); G.push_back(pav); }

using namespace std;

int main()
{
    vector<vector<int> > G;
    const int v0[] = {1,2,3,4}; push_array(G, v0);
    const int v1[] = {0,2}; push_array(G, v1);
    const int v2[] = {0,1,3}; push_array(G, v2);
    const int v3[] = {0,4,2}; push_array(G, v3);
    const int v4[] = {0,3}; push_array(G, v4);

    homsearch *h = new_homsearch(G, G, -1, true, true, -1);
    h->search(0);

    for (auto r: h->res_list) {
        for (auto fv: r)
            cout << fv << " ";
        cout << "\n";
    }
    assert(h->res_count == 6);

    delete h;

    return 0;
}
