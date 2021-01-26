// This code is part of the project "Theoretically Efficient Parallel Graph
// Algorithms Can Be Fast and Scalable", presented at Symposium on Parallelism
// in Algorithms and Architectures, 2018.
// Copyright (c) 2018 Laxman Dhulipala, Guy Blelloch, and Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all  copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// Usage:
// numactl -i all ./KCore -rounds 3 -s -m com-orkut.ungraph.txt_SJ
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -rounds : the number of times to run the algorithm
//     -fa : run the fetch-and-add implementation of k-core
//     -nb : the number of buckets to use in the bucketing implementation

// #define OUTPUT_EDGES
// #define OUTPUT_RESULT

#include "FullyDynamic.h"

// namespace gbbs {


// }  // namespace gbbs

void output_answer(gbbs::FullyDynamic fullyDynamic, const char fileName[], int nn = -1) {
	cout << "writing edges to output file" << endl;
	ofstream outputFile;
	outputFile.open(fileName);
	int n = nn;
	if(n == -1){
		n = fullyDynamic.maxNodeId();
	}
	for (int i = 1; i <= n; ++i)
		outputFile << i << " " << fullyDynamic.getApproxCoreVal(i) << endl;
	outputFile.close();
	cout << "finish writing" << endl;
}

int main(int argc, char **argv) {
	double epsilon = atof(argv[1]);
	double lambda = atof(argv[2]);
	double alpha = atof(argv[3]);
	char *fileName = argv[4];
	char *outFileName = argv[5];
	int upd = atoi(argv[6]); // new, 1 for insert, 2 for delete

    Update updtype;
	if(upd == 1){
		updtype = INS;
	}else{
		updtype = DEL;
	}

	gbbs::FullyDynamic fullyDynamic = gbbs::FullyDynamic(epsilon, lambda, alpha, fileName, updtype);

	time_t t = clock();
	fullyDynamic.run(INS);
	t = clock() - t;
	FILE *ofp = fopen("StatTimeMemoryIns.txt", "a");
	fprintf(ofp, "Round\t%.2f\t%f\t", epsilon, t / 1000.0);
	cerr << t << " ms." << endl;
	cerr << "Insertion Finished!!!" << endl;	
	fclose(ofp);

#ifdef OUTPUT_RESULT
	output_answer(fullyDynamic, outFileName);
#endif

	fullyDynamic.run(DEL);
	t = clock() - t;
	ofp = fopen("StatTimeMemoryDel.txt", "a");
	fprintf(ofp, "Round\t%.2f\t%f\t", epsilon, t / 1000.0);
	cerr << t << " ms." << endl;
	cerr << "Deletiontion Finished!!!" << endl;
	fclose(ofp);
	return 0;
}
