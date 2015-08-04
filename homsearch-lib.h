#include <bitset>
#include <vector>
#include <cassert>
#include <cstdint>

using namespace std;

#define MAX_GRAPH_SIZE 256

template< size_t size_lim >
class homsearch;

template< size_t size_lim >
class search_state {
  public:
    // Partial map, -1 for unmapped vertices
    vector<int> f;

    // Candidate targets for every vertex
    vector<bitset<size_lim> > candidates;
    vector<int> candidates_c;

    // Convenience pointer
    const homsearch<size_lim> *search;

  public:
    search_state(const homsearch<size_lim> *search_, const vector<int> *f_ = NULL):
      f(search_->G.size()), candidates(search_->G.size()), candidates_c(search_->G.size()), search(search_)
    {
        assert((f_ == NULL) || (f_->size() == G.size()));

        // Initialize full candidate lists
        for (unsigned int v = 0; v < search_->G.size(); v++) {
            for (unsigned int i = 0; i < search_->H.size(); i++)
                  candidates[v][i] = 1;
            candidates_c[v] = search_->H.size();
        }

        // Set partial map and limit candidates
        for (unsigned int i = 0; i < search->G.size(); i++) {
            if (f_ && ((*f_)[i] != -1)) {
                set_map(i, (*f_)[i]);
            } else {
                f[i] = -1;
            }
        }
    }

    // Limit the candidates of v to neighbors of n
    void inline limit_to_neighbors(int v, int n)
    {
	candidates[v] &= search->H_neighbors[n];
	candidates_c[v] = candidates[v].count();
    }

    // Set one vertex map and limit neighbor candidates
    void inline set_map(int v, int fv)
    {
	f[v] = fv;
	for (int n: search->G[v]) {
	    limit_to_neighbors(n, fv);
	}
    }

};

template< size_t size_lim >
class homsearch {
   public:
    // Graphs and neighbourhoods - external arrays
    const vector<vector<int> > G;
    vector <bitset<size_lim> > G_neighbors;

    const vector<vector<int> > H;
    vector<bitset<size_lim> > H_neighbors;

    // Results
    long long int res_limit;
    long long int res_count;
    vector<vector<int> > res_list;
    bool res_store;

    int max_depth;
    bool retract_mode;

   public:
    homsearch(const vector<vector<int> > G_, const vector<vector<int> > H_,
              long long int res_limit_, bool res_store_, bool retract_mode_, int max_depth_ = -1):
      G(G_), G_neighbors(G.size()), H(H_), H_neighbors(H.size()),
      res_count(0), res_limit(res_limit_), res_list(), res_store(res_store_),
      retract_mode(retract_mode_), max_depth(max_depth_)
    {   
        // Neighbor map in G
        for (int v = 0; v < G.size(); v++)
            for (auto i: G[v])
                G_neighbors[v][i] = 1;

        // Neighbor map in H
        for (int v = 0; v < H.size(); v++)
            for (auto i: H[v])
                H_neighbors[v][i] = 1;
    }

    homsearch(const homsearch &from):
      G(from.G), G_neighbors(from.G_neighbors), H(from.H), H_neighbors(from.H_neighbors),
      res_count(0), res_limit(from.res_limit), res_list(), res_store(from.res_store),
      max_depth(from.max_depth), retract_mode(from.retract_mode) { }


   public:
    void search(const search_state<size_lim> &s, int depth = 0);

   protected:
    void add_res(const search_state<size_lim> &s)
    {
        if (res_store && (res_count < res_limit))
            res_list.push_back(s.f);
        res_count ++;
    }
};


template< size_t size_lim >
void homsearch<size_lim>::search(const search_state<size_lim> &s, int depth)
{
    // Select branching vertex minimizing #candidates and
    int v = -1;
    int min_cand = H.size() + 1;
    int max_deg = -1;

    for (unsigned int i = 0; i < G.size(); i ++) {
	if ((s.f[i] == -1) && (s.candidates_c[i] <= min_cand)) {
	    if ((s.candidates_c[i] < min_cand) || (G[i].size() > max_deg)) {
		max_deg = G[i].size();
		min_cand = s.candidates_c[i];
		v = i;
	    }
	}
    }
    assert(min_cand > 0);

    // All vertices have been mapped
    if (v == -1) {
        add_res(s);
        return;
    }

    // Go over the candidates for v
    for (unsigned int fv = 0; (fv < H.size()); fv ++) {
        if (res_count >= res_limit) break;
	if (! s.candidates[v][fv]) continue;

	// Create subsearch
        search_state<size_lim> s2 = s;
	s2.set_map(v, fv);

	// Handle retract/core heuristics
	assert((s2.f[fv] == -1) || (s2.f[fv] == fv));
        if (s2.f[fv] == -1) {
            s2.set_map(fv, fv);
        }
        if (fv != v)
            for (unsigned int i = 0; i < G.size(); i ++)
                if (s2.f[i] == -1) {
                    s2.candidates[i][v] = 0;
                    s2.candidates_c[i] = s2.candidates[i].count();
                }

	// Run subsearch
        if ((max_depth >= 0) && (depth < max_depth))
            search(s2, depth + 1);
        else
            add_res(s2);
    }
}




/*

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
	if (res_list) {
	    vector<int> vec(G_size, -1);
	    for (int i = 0; i < G_size; i ++)
		vec[i] = f[i];
	    res_list->push_back(vec);
	}
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
*/

typedef search_state<32> search_state32;
typedef homsearch<32> homsearch32;

typedef search_state<256> search_state256;
typedef homsearch<256> homsearch256;

