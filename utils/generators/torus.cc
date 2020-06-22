// Usage: ./gen_torus [-w] <output_file>
// Flags:
//   optional:
//     -w <to generate a weighted torus>
//
// ex:
// > ./gen_torus 1000000000 3Dtorus_1B.adj
// > ./gen_torus -w 1000000000 3Dtorus_1B_wgh.adj

// Note: the output graph is stored in the uncompressed format. If you want to
// generate a massive torus graph and store it in the compressed format you
// should implement 3D-hilbert-ordering; currently the ids are not laid out
// optimally and will result in poor compression.

#include <math.h>
#include "../src/oldlib/benchIO.h"
#include "../src/oldlib/parse_command_line.h"
#include "gbbs/gbbs.h"
#include "gbbs/utils.h"

#include <algorithm>
#include <cstring>
#include <iostream>

#include "lib/random.h"

using namespace std;

namespace gbbs {

constexpr int max_weight = 32;

// rotate/flip a quadrant appropriately
void rot(long n, long* x, long* y, long rx, long ry) {
  if (ry == 0) {
    if (rx == 1) {
      *x = n - 1 - *x;
      *y = n - 1 - *y;
    }

    // Swap x and y
    long t = *x;
    *x = *y;
    *y = t;
  }
}

long xy2d(long n, long x, long y) {
  long rx, ry, s, d = 0;
  for (s = n / 2; s > 0; s /= 2) {
    rx = (x & s) > 0;
    ry = (y & s) > 0;
    d += s * s * ((3 * rx) ^ ry);
    rot(s, &x, &y, rx, ry);
  }
  return d;
}

uintE loc2d(uintE n, uintE i1, uintE i2) {
  i1 = (i1 + n) % n;
  i2 = (i2 + n) % n;
  return xy2d(n * n, i1, i2);
  //  return ((i1 + n) % n)*n + (i2 + n) % n;
}

// public static uint[] HilbertIndexTransposed(this uint[] hilbertAxes, int
// bits)
//{
//    var X = (uint[])hilbertAxes.Clone();
//    var n = hilbertAxes.Length; // n: Number of dimensions
//    uint M = 1U << (bits - 1), P, Q, t;
//    int i;
//    // Inverse undo
//    for (Q = M; Q > 1; Q >>= 1)
//    {
//        P = Q - 1;
//        for (i = 0; i < n; i++)
//            if ((X[i] & Q) != 0)
//                X[0] ^= P; // invert
//            else
//            {
//                t = (X[0] ^ X[i]) & P;
//                X[0] ^= t;
//                X[i] ^= t;
//            }
//    } // exchange
//    // Gray encode
//    for (i = 1; i < n; i++)
//        X[i] ^= X[i - 1];
//    t = 0;
//    for (Q = M; Q > 1; Q >>= 1)
//        if ((X[n - 1] & Q)!=0)
//            t ^= Q - 1;
//    for (i = 0; i < n; i++)
//        X[i] ^= t;
//
//    return X;
//}

uintE loc3d(uintE n, uintE i1, uintE i2, uintE i3) {
  return ((i1 + n) % n) * n * n + ((i2 + n) % n) * n + (i3 + n) % n;
}

template <class T>
void writeArrayToStream(ofstream& os, T* A, long n) {
  long BSIZE = 1000000;
  long offset = 0;
  while (offset < n) {
    // Generates a string for a sequence of size at most BSIZE
    // and then wrties it to the output stream
    std::cout << "Writing offset = " << offset << std::endl;
    _seq<char> S = benchIO::arrayToString(A + offset, min(BSIZE, n - offset));
    os.write(S.A, S.n);
    S.del();
    offset += BSIZE;
  }
}

void graph2DTorus(uintE n, char* fname) {
  uintE dn = round(pow((float)n, 1.0 / 2.0));
  uintE nn = dn * dn;
  size_t deg = 4;
  size_t nonZeros = deg * nn;
  std::cout << "nn = " << nn << " nz = " << nonZeros << " dn = " << dn << std::endl;
  uintE* edges = newA(uintE, nonZeros);
  parallel_for(size_t i = 0; i < dn; i++) {
    parallel_for(size_t j = 0; j < dn; j++) {
      uintE v = loc2d(dn, i, j);
      uintE nghs[4];
      nghs[0] = loc2d(dn, i, j + 1);
      nghs[1] = loc2d(dn, i, j - 1);
      nghs[2] = loc2d(dn, i + 1, j);
      nghs[3] = loc2d(dn, i - 1, j);
      std::sort(nghs, nghs + 4);
      size_t off = v * deg;
      for (size_t k = 0; k < deg; k++) {
        edges[off + k] = nghs[k];
      }
    }
  }
  uintT* degs = newA(uintT, nn);
  parallel_for(size_t i = 0; i < nn; i++) { degs[i] = i * 4; }

  ofstream file(fname, ios::out | ios::binary);
  if (!file.is_open()) {
    std::cout << "Unable to open file: " << fname << std::endl;
    exit(0);
  }
  file << "AdjacencyGraph" << std::endl;
  file << nn << std::endl;
  file << nonZeros << std::endl;
  writeArrayToStream(file, degs, nn);
  writeArrayToStream(file, edges, nonZeros);
  file.close();
  std::cout << "Wrote file." << std::endl;
  free(edges);
  free(degs);
}

void graph3DTorus(uintE n, char* fname) {
  uintE dn = round(pow((float)n, 1.0 / 3.0));
  uintE nn = dn * dn * dn;
  size_t deg = 6;
  size_t nonZeros = deg * nn;
  std::cout << "nn = " << nn << " nz = " << nonZeros << " dn = " << dn << std::endl;
  uintE* edges = newA(uintE, nonZeros);
  parallel_for(size_t i = 0; i < dn; i++) {
    parallel_for(size_t j = 0; j < dn; j++) {
      parallel_for(size_t k = 0; k < dn; k++) {
        uintE v = loc3d(dn, i, j, k);
        uintE nghs[6];
        nghs[0] = loc3d(dn, i, j, k + 1);
        nghs[1] = loc3d(dn, i, j, k - 1);
        nghs[2] = loc3d(dn, i, j + 1, k);
        nghs[3] = loc3d(dn, i, j - 1, k);
        nghs[4] = loc3d(dn, i + 1, j, k);
        nghs[5] = loc3d(dn, i - 1, j, k);
        std::sort(nghs, nghs + deg);
        size_t off = v * deg;
        for (size_t l = 0; l < deg; l++) {
          edges[off + l] = nghs[l];
        }
      }
    }
  }
  uintT* degs = newA(uintT, nn);
  parallel_for(size_t i = 0; i < nn; i++) { degs[i] = i * deg; }

  ofstream file(fname, ios::out | ios::binary);
  if (!file.is_open()) {
    std::cout << "Unable to open file: " << fname << std::endl;
    exit(0);
  }
  file << "AdjacencyGraph" << std::endl;
  file << nn << std::endl;
  file << nonZeros << std::endl;
  writeArrayToStream(file, degs, nn);
  writeArrayToStream(file, edges, nonZeros);
  file.close();
  std::cout << "Wrote file." << std::endl;
  free(edges);
  free(degs);
}

void graph3DTorusWgh(uintE n, char* fname) {
  uintE dn = round(pow((float)n, 1.0 / 3.0));
  uintE nn = dn * dn * dn;
  size_t deg = 6;
  size_t nonZeros = deg * nn;
  std::cout << "nn = " << nn << " nz = " << nonZeros << " dn = " << dn << std::endl;
  uintE* edges = newA(uintE, nonZeros);
  parallel_for(size_t i = 0; i < dn; i++) {
    parallel_for(size_t j = 0; j < dn; j++) {
      parallel_for(size_t k = 0; k < dn; k++) {
        uintE v = loc3d(dn, i, j, k);
        uintE nghs[6];
        nghs[0] = loc3d(dn, i, j, k + 1);
        nghs[1] = loc3d(dn, i, j, k - 1);
        nghs[2] = loc3d(dn, i, j + 1, k);
        nghs[3] = loc3d(dn, i, j - 1, k);
        nghs[4] = loc3d(dn, i + 1, j, k);
        nghs[5] = loc3d(dn, i - 1, j, k);
        std::sort(nghs, nghs + deg);
        size_t off = v * deg;
        for (size_t l = 0; l < deg; l++) {
          edges[off + l] = nghs[l];
        }
      }
    }
  }
  uintT* degs = newA(uintT, nn);
  parallel_for(size_t i = 0; i < nn; i++) { degs[i] = i * deg; }

  ofstream file(fname, ios::out | ios::binary);
  if (!file.is_open()) {
    std::cout << "Unable to open file: " << fname << std::endl;
    exit(0);
  }
  file << "AdjacencyGraph" << std::endl;
  file << nn << std::endl;
  file << nonZeros << std::endl;
  writeArrayToStream(file, degs, nn);
  writeArrayToStream(file, edges, nonZeros);
  file.close();
  std::cout << "Wrote file." << std::endl;
  free(edges);
  free(degs);
}

void graph3DTorus27(uintE n, char* fname) {
  uintE dn = round(pow((float)n, 1.0 / 3.0));
  uintE nn = dn * dn * dn;
  size_t deg = 26;
  size_t nonZeros = deg * nn;
  std::cout << "nn = " << nn << " nz = " << nonZeros << " dn = " << dn << std::endl;
  uintE* edges = newA(uintE, nonZeros);
  parallel_for(size_t i = 0; i < dn; i++) {
    parallel_for(size_t j = 0; j < dn; j++) {
      parallel_for(size_t k = 0; k < dn; k++) {
        size_t v = loc3d(dn, i, j, k);
        uintE nghs[26];
        uintE ct = 0;
        for (long l = -1; l < 2; l++) {
          for (long m = -1; m < 2; m++) {
            for (long n = -1; n < 2; n++) {
              uintE coord = loc3d(dn, i + l, j + m, k + n);
              if (coord != v) {
                nghs[ct++] = coord;
              }
            }
          }
        }
        assert(ct == 26);
        //        nghs[0] = loc3d(dn, i, j, k+1);
        //        nghs[1] = loc3d(dn, i, j, k-1);
        //        nghs[2] = loc3d(dn, i, j+1, k);
        //        nghs[3] = loc3d(dn, i, j-1, k);
        //        nghs[4] = loc3d(dn, i+1, j, k);
        //        nghs[5] = loc3d(dn, i-1, j, k);
        std::sort(nghs, nghs + deg);
        size_t off = v * deg;
        for (size_t l = 0; l < deg; l++) {
          edges[off + l] = nghs[l];
        }
      }
    }
  }
  uintT* degs = newA(uintT, nn);
  parallel_for(size_t i = 0; i < nn; i++) { degs[i] = i * deg; }

  ofstream file(fname, ios::out | ios::binary);
  if (!file.is_open()) {
    std::cout << "Unable to open file: " << fname << std::endl;
    exit(0);
  }
  file << "AdjacencyGraph" << std::endl;
  file << nn << std::endl;
  file << nonZeros << std::endl;
  std::cout << "writing m = " << nonZeros << std::endl;
  writeArrayToStream(file, degs, nn);
  writeArrayToStream(file, edges, nonZeros);
  file.close();
  std::cout << "Wrote file." << std::endl;
  free(edges);
  free(degs);
}

void BuildTorus(int argc, char* argv[]) {
  commandLine P(argc, argv, "[-w] n <outFile>");
  pair<int, char*> in = P.sizeAndFileName();
  long n = in.first;
  char* fname = in.second;
  bool weighted = P.getOptionValue("-w");
  std::cout << "Generating 3D torus" << std::endl;
  if (weighted) {
    graph3DTorusWgh(n, fname);
  } else {
    graph3DTorus(n, fname);
  }
}

}  // namespace gbbs

int main(int argc, char* argv[]) {
  BuildTorus(argc, argv);
}

