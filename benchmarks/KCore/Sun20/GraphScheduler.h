#ifndef __GRAPHSCHEDULER__
#define __GRAPHSCHEDULER__

#include <vector>
#include <fstream>
#include <cassert>
#include <string>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include "Hypergraph.h"
#include "gbbs/gbbs.h"
using namespace std;


enum Update {INS, DEL};

struct EdgeUpdate {
	Hyperedge e;
	int timestamp;
	Update updType;
};

template<class Graph>
inline vector<EdgeUpdate> shuffle_edges(Graph G, Update up){
  	vector<EdgeUpdate> updates_shuffled;
    size_t m = G.num_edges();
    auto perm = pbbs::random_permutation<gbbs::uintE>(m/2);
    pbbs::sequence<std::tuple<gbbs::uintE, gbbs::uintE, typename Graph::weight_type>> edge_list = G.edges();
    auto edge_list_dedup = pbbs::filter(edge_list, [&](const std::tuple<gbbs::uintE, gbbs::uintE, typename Graph::weight_type> & e){
      return std::get<0>(e) < std::get<1>(e);
    });
    edge_list.clear();
    // assert(edge_list_dedup.size() == m/2);
    for (size_t i = 0; i< edge_list_dedup.size(); ++i) {
		EdgeUpdate edgeUpdate;
		gbbs::uintE u = std::get<0>(edge_list_dedup[perm[i]]);
		gbbs::uintE v = std::get<1>(edge_list_dedup[perm[i]]);
		edgeUpdate.e.push_back(u);
		edgeUpdate.e.push_back(v);
		edgeUpdate.timestamp = i;
		edgeUpdate.updType = up;	
      	updates_shuffled.emplace_back(edgeUpdate);
    }
    edge_list_dedup.clear();
    std::cout << "shuffled and deduped" << std::endl;
    std::cout << updates_shuffled.size() << " deduped edges" << std::endl;
    return updates_shuffled;    
}

  template<class Graph>
  size_t get_max_deg(Graph& DG) {
    size_t max_deg = 0;
    gbbs::parallel_for(0, DG.n, [&] (size_t i) {
      size_t deg = DG.get_vertex(i).out_degree();
      pbbs::write_min(&max_deg, deg, std::greater<size_t>());
    });
    return max_deg;
  }

class GraphScheduler {
public:
	GraphScheduler(const char fileName[], Update t_updtype){
	updtype = t_updtype;
	// fin.open(fileName);
	G = gbbs::gbbs_io::read_unweighted_symmetric_graph(fileName, false);
	load();
	// fin.close();
	G.del();
	position = 0;
	}
	
	EdgeUpdate nextUpdate(){
		return updates[position++];
	}
	
	inline bool hasNext() {
		return position < updates.size();
	}
	unsigned numberOfNodes, maxDegree;
	Update updtype;
private:
	void load(){
		updates = shuffle_edges(G, updtype);

		// Only in C++11
		updates.shrink_to_fit();
		//edge_queue_no_time_.shrink_to_fit();
		numberOfNodes = G.num_vertices();
		maxDegree = get_max_deg(G);

		cerr << "Finished. " << updates.size() << " updates." << endl;
	}
	std::vector<EdgeUpdate> updates;
	// std::ifstream fin;
	gbbs::symmetric_graph<gbbs::symmetric_vertex, pbbslib::empty> G;
	unsigned position;
};
#endif
