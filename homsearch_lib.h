
#ifndef _HOMSEARCH_LIB_H_
#define _HOMSEARCH_LIB_H_

#include <bitset>
#include <vector>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <algorithm>

using namespace std;


/////////////////////////////////////
// Generic interface - virtual class

class homsearch {
   public:
    // Graphs
    const vector<vector<int> > G;
    const vector<vector<int> > H;

    // Results
    long long int res_limit;
    long long int res_count;
    vector<vector<int> > res_list;
    bool res_store;

    // Options
    int max_depth;
    bool retract_mode;

   public:
    homsearch(const vector<vector<int> > &G_, const vector<vector<int> > &H_,
              long long int res_limit_, bool res_store_, bool retract_mode_, int max_depth_):
      G(G_), H(H_), 
      res_limit(res_limit_), res_count(0), res_list(), res_store(res_store_),
      max_depth(max_depth_), retract_mode(retract_mode_) {}

    homsearch(const homsearch &from):
      G(from.G), H(from.H),
      res_limit(from.res_limit), res_count(0), res_list(), res_store(from.res_store),
      max_depth(from.max_depth), retract_mode(from.retract_mode) {}

    virtual ~homsearch() = default;

   public:

    virtual void search_vector(const vector<int> &f, int depth = 0) = 0;

    virtual void search(int depth = 0)
    {
        vector<int> f0(G.size(), -1);
        search_vector(f0, depth);
    }

};



////////////////////////////////////////////////
// Search state struct for particular set sizes

template< size_t size_lim >
class homsearch_impl;

template< size_t size_lim >
class homsearch_state {
  public:
    // Partial map, -1 for unmapped vertices
    vector<int> f;

    // Candidate targets for every vertex
    vector<bitset<size_lim> > candidates;

    // Convenience pointer
    const homsearch_impl<size_lim> *search;

  public:
    homsearch_state(const homsearch_impl<size_lim> *search_, const vector<int> *f_ = NULL):
      f(search_->G.size()), candidates(search_->G.size()), search(search_)
    {
        assert((f_ == NULL) || (f_->size() == search->G.size()));

        // Initialize full candidate lists
        for (unsigned int v = 0; v < search_->G.size(); v++) {
            for (unsigned int i = 0; i < search_->H.size(); i++)
                  candidates[v][i] = 1;
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

    // Set one vertex map, limit neighbor candidates and other heuristics
    // Returns success: if false, contradiction was found and mapping is not valid, state is broken
    bool inline set_map(int v, int fv)
    {
        assert(candidates[v][fv]);
        assert(f[v] == -1);
	f[v] = fv;

        // Find dist=1 vertices
        bitset<size_lim> N1G = search->G_neighbors[v];
        bitset<size_lim> N1H = search->H_neighbors[fv];

        // Limit dist=1 neighborhood candidates
        for (unsigned int n = 0; n < N1G.size(); n++)
            if ((f[n] == -1) && (N1G[n]))
                candidates[n] &= N1H;

        // Find dist=2 vertices
        bitset<size_lim> N2G;
        for (unsigned int n = 0; n < N1G.size(); n++)
            if (N1G[n])
                N2G |= search->G_neighbors[n];
        bitset<size_lim> N2H;
        for (unsigned int n = 0; n < N1H.size(); n++)
            if (N1H[n])
                N2H |= search->H_neighbors[n];
                
        // Limit dist=2 neighborhood candidates
        for (unsigned int n = 0; n < N2G.size(); n++)
            if ((f[n] == -1) && (N2G[n]))
                candidates[n] &= N2H;

        cout << "set_map(" << v << ", " << fv << ")\n    N1G=" << N1G << " N1H=" << N1H << "\n    N2G=" << N2G << " N2H=" << N2H << "\n";

        if (search->retract_mode) {

            // Retract/core heuristics

            // Target not mapped or fix-point (consistency of past candidates)
            assert((f[fv] == -1) || (f[fv] == fv));

            // Non-mapped target must be fix-point
            if (f[fv] == -1) {
                if (! candidates[fv][fv])
                    return false;
                if (! set_map(fv, fv))
                    return false;
            }

            // If this is not a fix-point, disable it as a target
            if (fv != v)
                for (unsigned int i = 0; i < candidates.size(); i ++)
                    if (f[i] == -1) {
                        candidates[i][v] = 0;
                    }
        } else {

            // Homomorphism heuristics

        }

        return true;
    }
};


//////////////////////////////////////////////////
// Implementations for particular fixed set sizes

template< size_t size_lim >
class homsearch_impl: public homsearch {

   public:
    vector <bitset<size_lim> > G_neighbors;
    vector<bitset<size_lim> > H_neighbors;

   public:
    homsearch_impl(const vector<vector<int> > &G_, const vector<vector<int> > &H_,
              long long int res_limit_, bool res_store_, bool retract_mode_, int max_depth_ = -1):
      homsearch(G_, H_, res_limit_, res_store_, retract_mode_, max_depth_),
      G_neighbors(G.size()), H_neighbors(H.size())
    {   
        // Neighbor map in G
        for (unsigned int v = 0; v < G.size(); v++)
            for (auto i: G[v])
                G_neighbors[v][i] = 1;

        // Neighbor map in H
        for (unsigned int v = 0; v < H.size(); v++)
            for (auto i: H[v])
                H_neighbors[v][i] = 1;
    }

    homsearch_impl(const homsearch_impl<size_lim> &from):
      homsearch(from),
      G_neighbors(from.G_neighbors), H_neighbors(from.H_neighbors) {}

    virtual ~homsearch_impl() = default;


   public:
    virtual void search_state(const homsearch_state<size_lim> &s, int depth = 0);

    virtual void search_vector(const vector<int> &f, int depth = 0)
    {
        homsearch_state<size_lim> s0(this, &f);
        search_state(s0, depth);
    }

   protected:
    void add_res(const homsearch_state<size_lim> &s)
    {
        if (res_store && ((res_count < res_limit) || (res_limit == -1)))
            res_list.push_back(s.f);
        res_count ++;
    }
};

template< size_t size_lim >
void homsearch_impl<size_lim>::search_state(const homsearch_state<size_lim> &s, int depth)
{
    // Select branching vertex minimizing #candidates and
    int v = -1;
    int min_cand = H.size() + 1;
    int max_deg = -1;

//    cout << "\nsearch d = " << depth << " / " << max_depth << "\n";
//    for (unsigned int i = 0; i < G.size(); i ++)
//        cout << i << " " << s.f[i] << " " << s.candidates[i] << "\n";

    for (unsigned int i = 0; i < G.size(); i ++) {
        int ccount = s.candidates[i].count();
	if ((s.f[i] == -1) && (ccount <= min_cand)) {
	    if ((ccount < min_cand) || ((int)G[i].size() > max_deg)) {
		max_deg = G[i].size();
		min_cand = ccount;
		v = i;
	    }
	}
    }

    // Some vertex has no candidates
    if (min_cand == 0)
        return;

    // All vertices have been mapped
    if (v == -1) {
        add_res(s);
        return;
    }

    // Go over the candidates for v
    for (unsigned int fv = 0; (fv < H.size()); fv ++) {
        if ((res_limit >= 0) && (res_count >= res_limit)) break;
	if (! s.candidates[v][fv]) continue;

	// Create subsearch
        homsearch_state<size_lim> s2(s);
	
        // Set map, check consistency
        if (! s2.set_map(v, fv))
            continue;

	// Run subsearch
        if ((max_depth >= 0) && (depth >= max_depth)) {
            add_res(s2);
        } else {
            search_state(s2, depth + 1);
        }
    }
}


///////////////////////////////////////////////////////////
// Helper to create the right instance of homsearch_impl<>

extern homsearch *new_homsearch(const vector<vector<int> > &G, const vector<vector<int> > &H,
              long long int res_limit, bool res_store, bool retract_mode, int max_depth=-1);

#endif // _HOMSEARCH_LIB_H_
