#ifndef __HYPERGRAPH__
#define __HYPERGRAPH__

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
using namespace std;

struct vectorHash {
	std::size_t operator()(std::vector<unsigned> const& vec) const {
		std::size_t seed = vec.size();
		for(auto& i: vec) {
			seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}
		return seed;
	}
};

typedef unsigned Node;
typedef std::vector<Node> Hyperedge;

class Hypergraph{
public:
	Hypergraph(){
		nNodes = nEdges = edgeIdCounter = 0;
	}
	
	unsigned nNodes;
	unsigned nEdges;
	unsigned edgeIdCounter;
	std::vector<Hyperedge> edgePool;
	std::unordered_multimap<Hyperedge, unsigned, vectorHash> edge2id;
	std::unordered_map<Node, std::unordered_set<unsigned> > eList;
	unsigned insertEdge(const Hyperedge& e){
		// Insert an edge to the hypergraph
		edge2id.insert(pair<const Hyperedge, unsigned>(e, edgeIdCounter));
		for (const Node u: e)
			eList[u].insert(edgeIdCounter);
		edgePool.push_back(e);
		// ++edgeIdCounter;
		++nEdges;
		return edgeIdCounter++;
	}

	void deleteEdge(const Hyperedge& e){
		// Delete an edge from the hypergraph
		auto iter = edge2id.find(e);
		for (const Node u: e)
			eList[u].erase(iter->second);
		edgePool[iter->second].clear();
		edge2id.erase(iter);
		--nEdges;
	}	
};
#endif
