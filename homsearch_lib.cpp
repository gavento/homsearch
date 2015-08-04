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


