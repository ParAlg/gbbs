import sys
import functools
import itertools
import time
import networkx as nx
from numpy import array




def cosimrank(G, u, v, alpha=0.9, 
             max_iter=100, tol=1.0e-6, weight='weight'):
  c = alpha
  alpha = 1
  if len(G) == 0:
    return {}

  if not G.is_directed():
    D = G.to_directed()
  else:
    D = G

  # Create a copy in (right) stochastic form
  W = nx.stochastic_graph(D, weight=weight)
  N = W.number_of_nodes()

  W_v = nx.stochastic_graph(D, weight=weight)
  N_v = W.number_of_nodes()

  x = dict.fromkeys(W, 0)
  x[u] = 1
  x_v = dict.fromkeys(W, 0)
  x_v[v] = 1

  p = dict.fromkeys(W, 1.0 / N)

  p_v = dict.fromkeys(W, 1.0 / N)

  dangling_weights = p
  dangling_weights_v = p_v
  dangling_nodes = [n for n in W if W.out_degree(n, weight=weight) == 0.0]
  dangling_nodes_v = [n for n in W if W.out_degree(n, weight=weight) == 0.0]

  sim = 0

  # power iteration: make up to max_iter iterations
  for iter in range(max_iter):
    xlast = x
    x = dict.fromkeys(xlast.keys(), 0)
    danglesum = alpha * sum(xlast[n] for n in dangling_nodes)
    for n in x:
      for nbr in W[n]:
        x[nbr] += alpha * xlast[n] * W[n][nbr][weight]
      x[n] += danglesum * dangling_weights.get(n, 0) + (1.0 - alpha) * p.get(n, 0)
    # check convergence, l1 norm
    err = sum([abs(x[n] - xlast[n]) for n in x])
    xlast_v = x_v
    x_v = dict.fromkeys(xlast_v.keys(), 0)
    danglesum_v = alpha * sum(xlast_v[n] for n in dangling_nodes_v)
    for n in x_v:
      # this matrix multiply looks odd because it is
      # doing a left multiply x^T=xlast^T*W
      for nbr in W_v[n]:
        x_v[nbr] += alpha * xlast_v[n] * W_v[n][nbr][weight]
      x_v[n] += danglesum_v * dangling_weights_v.get(n, 0) + (1.0 - alpha) * p_v.get(n, 0)
    # check convergence, l1 norm
    err_v = sum([abs(x_v[n] - xlast_v[n]) for n in x_v])
    sim += (c ** (iter+1)) * sum(x_v[key]*x.get(key, 0) for key in x_v)
    if err < N * tol and err_v < N_v * tol:
      return sim
  return sim





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

def ActualCoSimRank(G, u=0, v=1, importance_factor=0.85, max_iterations=100, tolerance=0.000001):
  t0 = time.time()
  similarity = cosimrank(G, u, v,alpha=importance_factor, max_iter=max_iterations, tol=tolerance)
  t1 = time.time()
  print("Time: ", t1-t0)
  print("Similarity: ", similarity)

def CoSimRank(G, src=0, ngh=1, importance_factor=0.85, max_iterations=100, tolerance=0.000001):
  print("Start SimRank")
  t0 = time.time()
  similarity = nx.simrank_similarity(G, source=src, target=ngh, importance_factor=importance_factor, max_iterations=max_iterations, tolerance=tolerance)
  t1 = time.time()
  print("Time: ", t1-t0)
  print("Similarity: ", similarity)

def CoSimRankNumpy(G, src=0, ngh=1, importance_factor=0.85, max_iterations=100, tolerance=0.0001):
  print("Start SimRankNumpy")
  t0 = time.time()
  similarity = nx.simrank_similarity_numpy(G, source=src, target=ngh, importance_factor=importance_factor, max_iterations=max_iterations, tolerance=tolerance)
  t1 = time.time()
  print("Time: ", t1-t0)
  print("Similarity: ", similarity)

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
  elif program_str == "CoSimRank":
    if argv_len > 4:
      CoSimRank(G, src=int(sys.argv[4]), ngh=int(sys.argv[5]))
    else:
      CoSimRank(G)
  elif program_str == "CoSimRankNumpy":
    if argv_len > 4:
      CoSimRank(G, src=int(sys.argv[4]), ngh=int(sys.argv[5]))
    else:
      CoSimRank(G)
  elif program_str == "ActualCoSimRank":
    if argv_len > 4:
      ActualCoSimRank(G, u=int(sys.argv[4]), v=int(sys.argv[5]))
    else:
      CoSimRank(G)

def main():
  argv_len = len(sys.argv)
  input_file = sys.argv[1] # First arg should be file name (ligra)
  symmetric = (sys.argv[2] == "s") # Second arg should be s or w/e for symmetric or not
  all_programs = ["ActualCoSimRank","CoSimRank","CoSimRankNumpy","BFS", "MaximalMatching", "KCore", "PageRank", "CliqueCounting", "TriangleCounting", "GeneralWeightSSSP", "GraphColoring"]
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