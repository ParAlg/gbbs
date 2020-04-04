// TODO(tomtseng): update this file comment
// Converts a SNAP graph (http://snap.stanford.edu/data/index.html) to
// Ligra adjacency graph format. To symmetrize the graph, pass the "-s"
// flag. For undirected graphs on SNAP, the "-s" flag must be passed
// since each edge appears in only one direction
//
//
// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#include <exception>

#include "ligra/parse_command_line.h"

int main(int argc, char* argv[]) {
  // TODO(tomtseng) write a more detailed help string
  const std::string kCommandLineHelpString{
    "[-s] [-w] -i <input SNAP file> -o <output AdjacencyGraph file>"};
  const std::string kDefaultOptionValue{""};
  const std::string kInputFlag{"-i"};
  const std::string kOutputFlag{"-o"};

  const commandLine parameters{argc, argv, kCommandLineHelpString};
  const std::string input_file{
    parameters.getOptionValue(kInputFlag, kDefaultOptionValue)};
  const std::string output_file{
    parameters.getOptionValue(kOutputFlag, kDefaultOptionValue)};
  const bool symmetrize{parameters.getOption("-s")};
  const bool weighted{parameters.getOption("-w")};

  if (input_file == kDefaultOptionValue || output_file == kDefaultOptionValue) {
    std::cerr << "ERROR: Please specify the input SNAP file with the '" <<
      kInputFlag << "' flag and the output file with the '" << kOutputFlag << 
      "' flag.\n";
    std::terminate();
  }

  if (weighted) {
    // TODO(tomtseng)
    (void)symmetrize;
    // wghEdgeArray<uintT> G = readWghSNAP<uintT>(iFile);
    // writeWghGraphToFile<uintT>(wghGraphFromWghEdges(G,sym),oFile);
  } else {
    // TODO(tomtseng)
    (void)symmetrize;
    // edgeArray<uintT> G = readSNAP<uintT>(iFile);
    // writeGraphToFile<uintT>(graphFromEdges(G,sym),oFile);
  }

  return 0;
}
