
#include "homsearch-lib.h"

template< std::size_t G_size_lim, std::size_t H_size_lim >
void homsearch<G_size_lim, H_size_lim>::search()
{
    int v = -1;

    // Select branching vertex
    int min_cand = H_size + 1;
    int max_deg = -1;

    for (int i = 0; i < H_size; i ++) {
	if ((f[i] != -1) && (candidates_c[i] <= min_cand)) {
	    if ((candidates_c[i] < min_cand) || (G_degs[i] > max_deg)) {
		max_deg = G_degs[i];
		min_cand = candidates_c[i];
		v = i;
	    }
	}
    }

    // All vertices mapped
    if (v == -1) {
	res_count ++;
	if (res_list)
	    res_list->append(f);
	return;
    }
	
	
    // Go over the candidates
    for (int fv = 0; (fv < H_size) && (res_limit > res_count); fv ++) {
	if (! candidates[v][fv]) continue;

	// Subsearch
	homsearch h2 = *this;
	h2.set_map(v, fv);

	// Handle retracting heuristics
	if (retract) {
	    assert((f[fv] == -1) || (f[fv] == fv));
	    if (f[fv] == -1) {
		h2.set_map(fv, fv);
	    }
	    if (fv != v)
		for (int i = 0; i < H_size; i ++)
		    if (f[i] == -1) {
			if (candidates[i][v])
			    candidates_c[i] --;
			candidates[i][v] = 0;
		    }
	}

	// Run subsearch
	h2.search();
	res_count = h2.res_count;
    }
}
