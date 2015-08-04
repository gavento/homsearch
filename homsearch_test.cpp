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

    homsearch *h = new_homsearch(G, G, 100, true, true, -1);
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
