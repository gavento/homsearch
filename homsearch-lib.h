#include <bitset>
#include <vector>
#include <cassert>
#include <cstdint>

using namespace std;

template< std::size_t G_size_lim, std::size_t H_size_lim >
class homsearch {
  public:
    const int G_size, H_size;
    // Graphs and neighbourhoods - external arrays
    const int *G_degs;
    const int **G_neighborlist;
    const bitset<H_size_lim> *H_neighbors;
    // Partial map, -1 for unmapped vertices
    int f[G_size_lim];
    // Candidate targets for every vertex
    bitset<H_size_lim> candidates[G_size_lim];
    int candidates_c[G_size_lim];
    // Results
    int64_t res_limit;
    int64_t res_count;
    vector<int[G_size_lim]> *res_list;
    // Looking for retracts only
    bool retract;

  public:
    void search();

    // Limit the candidates of v to neighbors of n
    void inline limit_candidates(int v, int n)
    {
	candidates[v] &= H_neighbors[n];
	candidates_c[v] = candidates[v].count();
    }

    void inline set_map(int v, int fv)
    {
	f[v] = fv;
	for (int i = 0; i < G_degs[v]; i ++){
	    limit_candidates(G_neighborlist[v][i], fv);
	}
    }
};

