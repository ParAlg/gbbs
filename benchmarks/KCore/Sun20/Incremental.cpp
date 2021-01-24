/*
Our insertion-only algorithm.
The code includes various functional blocks. Uncomment respective parts to examine different aspects of the execution.

[To compile]
g++ -std=c++11 -O3 Incremental.cpp GraphScheduler.cpp Hypergraph.cpp HypergraphCoreDecomp.cpp -o Incremental

[To run]
Incremental epsilon lambda filename

[Format of input]
The file should contain an update in each line.
Hyperedge insertion: + [node IDs] timestamp
So be careful that the last number of a line beginning with "+" is not an endpoint of the inserted hyperedge.
The timestamps are actually useless (for our purpose). You will need them only when you want to trace the insertion time of hyperedges.
*/

#include <cstdio>
#include <cmath>
#include <ctime>
#include <cassert>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include "Hypergraph.hpp"
#include "GraphScheduler.hpp"
#include "HypergraphCoreDecomp.hpp"
using namespace std;

class Incremental {
public:
	Incremental(double epsilon, double lambda, char fileName[]): epsilon(epsilon), lambda(lambda), scheduler(fileName) {
		initialize();
	}
	void run() {
		FILE *ofpVal = fopen("StatIncrementalCoreValue.txt", "w");
		FILE *ofpTime = fopen("StatIncrementalTime.txt", "w");
		FILE *ofpDetail = fopen("StatIncrementalDetail.txt", "w");
		int cnt = 0;
		time_t t0 = clock(), totalTime = 0;
		while (scheduler.hasNext()) {
			EdgeUpdate edgeUpdate = scheduler.nextUpdate();
			assert(edgeUpdate.updType == INS);
			insertEdge(edgeUpdate.e);
			++cnt;
			if (cnt % 100000 == 0) {
				fprintf(ofpTime, "%d\n", clock() - t0);
				fprintf(stderr, "%d...\t", cnt);
				fprintf(ofpVal, "%d\n", cnt);
				/*
				// Run the static exact algorithm
				HypergraphCoreDecomp hcd(h);
				hcd.solve();
				unordered_map<int, int> coreCount, coreCountApprox;
				for (auto& p: hcd.c) {
					const Node u = p.first;
					if (hcd.c[u] > 0) {
						if (!coreCount.count(hcd.c[u]))
							coreCount[hcd.c[u]] = 1;
						else
							++coreCount[hcd.c[u]];
						if (!coreCountApprox.count(b[tau][u]))
							coreCountApprox[b[tau][u]] = 1;
						else
							++coreCountApprox[b[tau][u]];
					}
				}
				for (auto& p: coreCount)
					fprintf(ofpVal, "%d %d\n", p.first, p.second);
				fprintf(ofpVal, "-1\n");
				for (auto& p: coreCountApprox)
					fprintf(ofpVal, "%d %d\n", p.first, p.second);
				fprintf(ofpVal, "-1\n");*/

				/*vector<double> maxErr(tau + 1, 0), avgErr(tau + 1, 0);
				int cntNonZero = 0;
				for (auto& p: hcd.c) {
					const Node u = p.first;
					if (hcd.c[u] > 0) {
						++cntNonZero;
						for (int t = 1; t <= tau; ++t) {
							double err = max(((double)hcd.c[u]) / b[t][u], ((double)b[t][u]) / hcd.c[u]);
							if (maxErr[t] < err)
								maxErr[t] = err;
							avgErr[t] += err;
						}
					}
				}
				for (int t = 1; t <= tau; ++t)
					avgErr[t] /= cntNonZero;
				for (int t = 1; t <= tau; ++t)
					fprintf(ofpDetail, "%.9f%c", maxErr[t], t == tau ? '\n' : ' ');
				for (int t = 1; t <= tau; ++t)
					fprintf(ofpDetail, "%.9f%c", avgErr[t], t == tau ? '\n' : ' ');
				*/
				t0 = clock();
			}
		}
	//	fclose(ofp);

		/*
		// Block: compare between different parameters
		FILE *ofpComp = fopen("StatIncremental.txt", "w");
		fprintf(ofpComp, "lambda = %f\n", lambda);
		fprintf(ofpComp, "total time = %d ms\n", clock() - t0 - totalExcludedTime);
		fprintf(ofpComp, "total number of updates = %d\n", cnt);
		double maxErr = 0, avgErr = 0;
		int cntNonZero = 0;
		HypergraphCoreDecomp hcd(h);
		hcd.solve();
		for (auto& p: hcd.c) {
			const Node u = p.first;
			assert((hcd.c[u] == 0 && b[tau][u] == 0) || (hcd.c[u] > 0 && b[tau][u] > 0));
			if (hcd.c[u] > 0) {
				++cntNonZero;
				double err = max((double)hcd.c[u] / b[tau][u], (double)b[tau][u] / hcd.c[u]);
				maxErr = max(maxErr, err);
				avgErr += err;
			}
		}
		avgErr /= cntNonZero;
		fprintf(ofpComp, "max error = %.9f\n", maxErr);
		fprintf(ofpComp, "avg error = %.9f\n", avgErr);
		fclose(ofpComp);
		*/

		/*
		// Block: show that a large tau is unnecessary
		vector<double> maxErr(tau + 1, 0), avgErr(tau + 1, 0);
		int cntNonZero = 0;
		HypergraphCoreDecomp hcd(h);
		hcd.solve();
		for (auto& p: hcd.c) {
			const Node u = p.first;
			assert((hcd.c[u] == 0 && b[tau][u] == 0) || (hcd.c[u] > 0 && b[tau][u] > 0));
			if (hcd.c[u] > 0) {
				++cntNonZero;
				for (int t = 1; t <= tau; ++t) {
					if (maxErr[t] < max((double)hcd.c[u] / b[t][u], (double)b[t][u] / hcd.c[u]))
						maxErr[t] = max((double)hcd.c[u] / b[t][u], (double)b[t][u] / hcd.c[u]);
					avgErr[t] += max((double)hcd.c[u] / b[t][u], (double)b[t][u] / hcd.c[u]);
				}
			}
		}
		for (int t = 1; t <= tau; ++t)
			avgErr[t] /= cntNonZero;

		FILE *ofpLargeTau = fopen("LargeTauIsUnnecessaryIncremental.txt", "w");
		fprintf(ofpLargeTau, "%d\n", tau);
		for (int t = 1; t <= tau; ++t)
			fprintf(ofpLargeTau, "%.9f %.9f\n", maxErr[t], avgErr[t]);
		fclose(ofpLargeTau);*/
	}
	unsigned getApproxCoreVal(Node u) {
		return b[tau][u];
	}
	Hypergraph h;
private:
	double epsilon, lambda;
	GraphScheduler scheduler;

	int tau;
	vector<unsigned> succ, pred;
	vector<unordered_map<Node, unsigned>> b, sigma;
	void initialize() {
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
			++sigma[1][u], bad.insert(u);
		for (unsigned t = 1; t <= tau; ++t) {
			bad2.clear();
			if (t < tau) {
				unsigned b_e = INT_MAX;
				for (const Node u: e)
					b_e = min(b_e, b[t][u]);
				for (const Node u: e) {
					if (b_e >= succ[b[t + 1][u]])
						++sigma[t + 1][u], bad2.insert(u);
				}
			}
			for (const Node u: bad)
				if (sigma[t][u] >= succ[b[t][u]])
					promote(t, u, bad2);
			swap(bad, bad2);
		}
	}
	void promote(const unsigned t, const Node u, unordered_set<Node>& bad2) {
	//	cerr << "Promote " << t << ' ' << u << ' ' << newEdgeId << endl;
		unsigned old_b_t_u = b[t][u];
		b[t][u] = succ[b[t][u]];
		updateSigma(t, u);
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
				if (old_b_e < succ[b[t + 1][v]]
							&& succ[b[t + 1][v]] <= new_b_e)
					++sigma[t + 1][v], bad2.insert(v);
			}
		}
	//	cerr << "Promote finished." << endl;
	}
	void updateSigma(unsigned t, Node u) {
		sigma[t][u] = 0;
		for (const unsigned eId: h.eList[u]) {
			const Hyperedge& e = h.edgePool[eId];
			unsigned b_e = INT_MAX;
			if (t > 1)
				for (const Node v: e)
					b_e = min(b_e, b[t - 1][v]);
			if (b_e >= succ[b[t][u]]) ++sigma[t][u];
		}
	}
};

int main(int argc, char **argv) {
	double epsilon = atof(argv[1]);
	double lambda = atof(argv[2]);
	char *fileName = argv[3];
	Incremental incremental(epsilon, lambda, fileName);
	time_t t = clock();
	incremental.run();
	cerr << (clock() - t) << " ms." << endl;
//	Node u;
//	while (cin >> u)
//		cout << incremental.getApproxCoreVal(u) << endl;
	return 0;
}
