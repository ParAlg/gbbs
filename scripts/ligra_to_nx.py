import sys
import functools
import itertools
import time
import networkx as nx

def read_ligra_symmetric_graph(input_file):
  G = nx.Graph()
  with open(input_file, 'r') as reader:
    reader.readline() # Discard AdjacencyGraph
    n = int(reader.readline())
    m = int(reader.readline())
    offsets = [None] * (n+1)
    offsets[n] = m
    edges = [(0, 0)] * m
    G.add_nodes_from(range(n))
    for i in range(n):
      offsets[i] = int(reader.readline())
    offsets_idx = 0
    for i in range(m):
      if offsets[offsets_idx + 1] == i:
        offsets_idx+=1
      edges[i] = (offsets_idx, int(reader.readline()))
    G.add_edges_from(edges)
  return G

def read_ligra_directed_graph(input_file):
  G = nx.DiGraph()
  with open(input_file, 'r') as reader:
    reader.readline() # Discard AdjacencyGraph
    n = int(reader.readline())
    m = int(reader.readline())
    offsets = [None] * (n+1)
    offsets[n] = m
    edges = [(0, 0)] * m
    G.add_nodes_from(range(n))
    for i in range(n):
      offsets[i] = int(reader.readline())
    offsets_idx = 0
    for i in range(m):
      if offsets[offsets_idx + 1] == i:
        offsets_idx+=1
      edges[i] = (offsets_idx, int(reader.readline()))
    G.add_edges_from(edges)
  return G

def TriangleCounting(G):
  print("Start TriangleCounting")
  t0 = time.time()
  tri_dict = nx.triangles(G)
  t1 = time.time()
  print("Time: ", t1-t0)
  total = 0
  for key in tri_dict:
    total += tri_dict[key]
  total = total / 3
  print("# Triangles: ",total)

def CliqueCounting(G):
  # TODO: only able to find all maximal cliques, not for specific k
  print("Start CliqueCounting")
  t0 = time.time()
  clique_iterator = nx.algorithms.clique.find_cliques(G)
  t1 = time.time()
  print("Time: ", t1-t0)

def PageRank(G, max_iter=100, tol=1e-06):
  print("Start PageRank")
  t0 = time.time()
  rank_dict = nx.pagerank(G, max_iter=max_iter, tol=tol)
  t1 = time.time()
  print("Time: ", t1-t0)
  max_pr = 0
  for key in rank_dict:
    max_pr = rank_dict[key] if rank_dict[key] > max_pr else max_pr
  print("Max PageRank: ",max_pr)

def KCore(G):
  print("Start KCore")
  t0 = time.time()
  kcore_subgraph = nx.algorithms.core.k_core(G)
  t1 = time.time()
  print("Time: ", t1-t0)
  print("Max core size: ",len(kcore_subgraph))

def MaximalMatching(G):
  print("Start MaximalMatching")
  t0 = time.time()
  # Set of edges
  matching_set = nx.algorithms.matching.maximal_matching(G)
  t1 = time.time()
  print("Time: ", t1-t0)

def BFS(G, src=0):
  print("Start BFS")
  t0 = time.time()
  edges_generator = nx.algorithms.traversal.breadth_first_search.bfs_edges(G, src)
  t1 = time.time()
  print("Time: ", t1-t0)

def GeneralWeightSSSP(G, src=0):
  print("Start GeneralWeightSSSP")
  t0 = time.time()
  # Dictionary keyed by targets with a list of nodes on shortest path from source
  paths_dict = nx.algorithms.shortest_paths.generic.shortest_path(G, source=src, method='bellman-ford')
  t1 = time.time()
  print("Time: ", t1-t0)

# strategy = 'random_sequential', 'smallest_last', 'independent_set', 
# 'connected_sequential_bfs', 'connected_sequential_dfs', 
# 'saturation_largest_first'
def GraphColoring(G, strategy='largest_first'):
  print("Start GraphColoring")
  t0 = time.time()
  coloring_dict = nx.algorithms.coloring.greedy_color(G, strategy=strategy)
  t1 = time.time()
  print("Time: ", t1-t0)

def program_parser(G, program_str):
  argv_len = len(sys.argv)
  if program_str == "BFS":
    if argv_len > 4:
      BFS(G, int(sys.argv[4]))
    else:
      BFS(G)
  elif program_str == "MaximalMatching":
    MaximalMatching(G)
  elif program_str == "KCore":
    KCore(G)
  elif program_str == "PageRank":
    if argv_len > 4:
      PageRank(G, int(sys.argv[4]), float(sys.argv[5]))
    else:
      PageRank(G)
  elif program_str == "CliqueCounting":
    CliqueCounting(G)
  elif program_str == "TriangleCounting":
    TriangleCounting(G)
  elif program_str == "GeneralWeightSSSP":
    if argv_len > 4:
      GeneralWeightSSSP(G, int(sys.argv[4]))
    else:
      GeneralWeightSSSP(G)
  elif program_str == "GraphColoring":
    if argv_len > 4:
      GraphColoring(G, sys.argv[4])
    else:
      GraphColoring(G)

def main():
  argv_len = len(sys.argv)
  input_file = sys.argv[1] # First arg should be file name (ligra)
  symmetric = (sys.argv[2] == "s") # Second arg should be s or w/e for symmetric or not
  all_programs = ["BFS", "MaximalMatching", "KCore", "PageRank", "CliqueCounting", "TriangleCounting", "GeneralWeightSSSP", "GraphColoring"]
  # use read edge list for snap format (TODO)
  G = read_ligra_symmetric_graph(input_file) if symmetric else read_ligra_directed_graph(input_file)
  program_str = sys.argv[3] # Third arg should be name of benchmark
  program_parser(G, program_str)


if __name__ == "__main__":
  main()


# ApproximateDensestSubgraph
# DegeneracyOrder
# Spanner
# ApproximateSetCover
# MinimumSpanningForest
# SpanningForest
# StronglyConnectedComponents
# Biconnectivity
# IntegralWeightSSSP
# PositiveWeightSSSP
# CoSimRank
# KTruss
# SCAN
# Connectivity
# LowDiameterDecomposition
# SSBetweenessCentrality
# CycleCounting
# MaximalIndependentSet
# SSWidestPath