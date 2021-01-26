/*
Our fully dynamic algorithm.
The code includes various functional blocks. Uncomment respective parts to examine different aspects of the execution.

[To compile]
g++ -std=c++11 -O3 FullyDynamic.cpp GraphScheduler.cpp Hypergraph.cpp HypergraphCoreDecomp.cpp -o FullyDynamic -lpsapi

[To run]
FullyDynamic epsilon lambda alpha filename

[Format of input]
The file should contain an update in each line.
Hyperedge insertion: + [node IDs] timestamp
Hyperedge deletion: - [node IDs]
So be careful that the last number of a line beginning with "+" is not an endpoint of the inserted hyperedge.
The timestamps are actually useless (for our purpose). You will need them only when you want to trace the insertion time of hyperedges.
When deleting a hyperedge, make sure that exactly the same sequence appeared in an insertion update before.
Therefore, it is recommended that in each line, the node IDs are sorted in some particular order, e.g., increasing order.

[Remark]
This program is developed on Windows. "-lpsapi" is related to reporting memory usage. You may have to remove it and "OutputMemory.cpp" (which includes a function outputMemory()) and reimplement this part when you want to compile and run the code on a different operating system.
*/

#include <cstdio>
#include <cmath>
#include <ctime>
#include <cassert>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include "gbbs/gbbs.h"
#include "Hypergraph.h"
#include "GraphScheduler.h"
#include "HypergraphCoreDecomp.h"
// #include "OutputMemory.cpp"
using namespace std;
namespace gbbs {




class FullyDynamic {
public:
	Hypergraph h;

	FullyDynamic(double epsilon, double lambda, double alpha, char fileName[], Update updtype):
		epsilon(epsilon), lambda(lambda), alpha(alpha), scheduler(fileName, updtype) {
		initialize();
	}

	unsigned maxNodeId(){
		return scheduler.maxNodeId;
	}
	// void run() {
	// 	FILE *ofpTime = fopen("StatFullyDynamicTime.txt", "w");
	// 	int cnt = 0;
	// 	time_t t0 = clock(), totalTime = 0;
	// 	int lastUpdTimestamp;
	// 	while (scheduler.hasNext()) {
	// 		EdgeUpdate edgeUpdate = scheduler.nextUpdate();
	// 		if (edgeUpdate.updType == INS)
	// 			insertEdge(edgeUpdate.e);
	// 		else
	// 			deleteEdge(edgeUpdate.e);
	// 		++cnt;

	// 		if (cnt % 100000 == 0) {
	// 			cerr << cnt << "...\t";
	// 			// Block: print all b[t][u] and c[u] every 100000 updates
	// 			time_t t1 = clock();
	// 			fprintf(ofpTime, "%d\n", t1 - t0);
	// 			t0 = clock();
	// 		}
	// 	} // end of while loop
	// 	fclose(ofpTime);
	// } //end of run function

	// if insert, inserting from 0 to end
	// if deleting, deleting from end to 0
	void run(Update updType) {
		int cnt = 0;
		time_t t0 = clock(), totalTime = 0;
		
		if(updType == INS){
			while (scheduler.hasNext()) {
				EdgeUpdate edgeUpdate = scheduler.nextUpdate();
				insertEdge(edgeUpdate.e);
				++cnt;

				if (cnt % 100000 == 0) {
					cerr << cnt << "...\t";
				}
			} // end of while loop
		}else{
			while (scheduler.hasPrev()) {
				EdgeUpdate edgeUpdate = scheduler.prevUpdate();
				deleteEdge(edgeUpdate.e);
				++cnt;

				if (cnt % 100000 == 0) {
					cerr << cnt << "...\t";
				}
			} // end of while loop
		}
	} //end of run function


	unsigned getApproxCoreVal(Node u) {
		return b[tau][u];
	}
	// void debug() {
	// 	for (int i = 1; i <= 100; ++i)
	// 		cerr << getApproxCoreVal(i) << ' ';
	// 	cerr << endl;
	// }
	// Hypergraph h;
private:
	double epsilon, lambda, alpha;
	GraphScheduler scheduler;

	int tau;
	vector<unsigned> succ, pred;
	vector<unordered_map<Node, unsigned>> b, sigma, rho;
	void initialize() {
		// alpha should be r * (1 + 3 * epsilon) where r is maximum edge cardinality.
		// But we can slightly reduce it when r is too large.
		// alpha is decided by input.

		// alpha = 2 * (1 + 3 * epsilon);
		tau = ceil(0.15 * log(scheduler.numberOfNodes) / log(1.0 + epsilon));

		// Build succ and pred
		vector<unsigned> Lambda(1, 0);
		int i = 0;
		while (Lambda[i] <= scheduler.maxDegree) {
			Lambda.push_back(max((unsigned)(Lambda[i] * (1.0 + lambda)), Lambda[i] + 1));
			++i;
		}
		succ = pred = vector<unsigned>(Lambda[i] + 1);
		for (int j = 0; j < i; ++j)
			succ[Lambda[j]] = Lambda[j + 1], pred[Lambda[j + 1]] = Lambda[j];
		b.resize(tau + 1);
		sigma.resize(tau + 1);
		rho.resize(tau + 1);
		cerr << "Lambda: ";
		for (auto& lambda: Lambda) cerr << lambda << ' ';
		cerr << endl;
		cerr << "# of nodes = " << scheduler.numberOfNodes << endl;
		cerr << "tau = " << tau << endl;
		cerr << "max degree = " << scheduler.maxDegree << endl;
	}
	void insertEdge(const Hyperedge& e) {
	//	cerr << "Insert an edge" << endl;
		h.insertEdge(e);
		unordered_set<Node> bad, bad2;
		for (const Node u: e)
			++sigma[1][u], ++rho[1][u], bad.insert(u);
		for (unsigned t = 1; t <= tau; ++t) {
			bad2.clear();
			if (t < tau) {
				unsigned b_e = INT_MAX;
				for (const Node u: e)
					b_e = min(b_e, b[t][u]);
				for (const Node u: e) {
					if (b_e >= succ[b[t + 1][u]])
						++sigma[t + 1][u], bad2.insert(u);
					if (b_e >= b[t + 1][u])
						++rho[t + 1][u];
				}
			}
			for (const Node u: bad)
				if (sigma[t][u] >= (unsigned)(alpha * succ[b[t][u]]))
					promote(t, u, bad2);
			swap(bad, bad2);
		}
	}
	void promote(const unsigned t, const Node u, unordered_set<Node>& bad2) {
	//	cerr << "Promote " << t << ' ' << u << ' ' << newEdgeId << endl;
		unsigned old_b_t_u = b[t][u];
		b[t][u] = succ[b[t][u]];
		updateSigmaAndRho(t, u);
	//	cerr << "b[" << t << "][" << u << "] = " << b[t][u] << ", sigma[" << t << "][" << u << "] = " << sigma[t][u] << endl;
		if (t == tau) return;
		for (const unsigned eId: h.eList[u]) {
	//		cerr << "eId = " << eId << endl;
			const Hyperedge& e = h.edgePool[eId];
			unsigned old_b_e = INT_MAX;
			for (const Node v: e)
				if (v != u)
					old_b_e = min(old_b_e, b[t][v]);
			unsigned new_b_e = min(old_b_e, b[t][u]);
			old_b_e = min(old_b_e, old_b_t_u);
			if (new_b_e == old_b_e) continue;
			for (const Node v: e) {
				if (old_b_e < succ[b[t + 1][v]] && succ[b[t + 1][v]] <= new_b_e)
					++sigma[t + 1][v], bad2.insert(v);
				if (old_b_e < b[t + 1][v] && b[t + 1][v] <= new_b_e)
					++rho[t + 1][v];
			}
		}
	//	cerr << "Promote finished." << endl;
	}
	void deleteEdge(const Hyperedge& e) {
	//	cerr << "Delete an edge" << endl;
		h.deleteEdge(e);
		unordered_set<Node> bad, bad2;
		for (const Node u: e)
			--sigma[1][u], --rho[1][u], bad.insert(u);
		for (unsigned t = 1; t <= tau; ++t) {
			bad2.clear();
			if (t < tau) {
				unsigned b_e = INT_MAX;
				for (const Node u: e)
					b_e = min(b_e, b[t][u]);
				for (const Node u: e) {
					if (b_e >= succ[b[t + 1][u]])
						--sigma[t + 1][u];
					if (b_e >= b[t + 1][u])
						--rho[t + 1][u], bad2.insert(u);
				}
			}
			for (const Node u: bad)
				if (rho[t][u] < b[t][u])
					demote(t, u, bad2);
			swap(bad, bad2);
		}
	}
	void demote(const unsigned t, const Node u, unordered_set<Node>& bad2) {
	//	cerr << "Demote " << t << ' ' << u << endl;
		unsigned old_b_t_u = b[t][u];
		b[t][u] = pred[b[t][u]];
		updateSigmaAndRho(t, u);
		if (t == tau) return;
		for (const unsigned eId: h.eList[u]) {
			const Hyperedge& e = h.edgePool[eId];
			unsigned old_b_e = INT_MAX;
			for (const Node v: e)
				if (v != u)
					old_b_e = min(old_b_e, b[t][v]);
			unsigned new_b_e = min(old_b_e, b[t][u]);
			old_b_e = min(old_b_e, old_b_t_u);
			if (new_b_e == old_b_e) continue;
			for (const Node v: e) {
				if (new_b_e < succ[b[t + 1][v]] && succ[b[t + 1][v]] <= old_b_e)
					--sigma[t + 1][v];
				if (new_b_e < b[t + 1][v] && b[t + 1][v] <= old_b_e)
					--rho[t + 1][v], bad2.insert(v);
			}
		}
	}
	void updateSigmaAndRho(unsigned t, Node u) {
		sigma[t][u] = rho[t][u] = 0;
		for (const unsigned eId: h.eList[u]) {
			const Hyperedge& e = h.edgePool[eId];
			unsigned b_e = INT_MAX;
			if (t > 1)
				for (const Node v: e)
					b_e = min(b_e, b[t - 1][v]);
			if (b_e >= succ[b[t][u]]) ++sigma[t][u];
			if (b_e >= b[t][u]) ++rho[t][u];
		}
	}
};

}  // namespace gbbs
