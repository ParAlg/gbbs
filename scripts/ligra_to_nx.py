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

def main():
  argv_len = len(sys.argv)
  input_file = sys.argv[1] # First arg should be file name (ligra)
  symmetric = (sys.argv[2] == "s") # Second arg should be s or w/e for symmetric or not
  # use read edge list for snap format (TODO)
  G = read_ligra_symmetric_graph(input_file) if symmetric else read_ligra_directed_graph(input_file)
  print("Start triangle count")
  t0 = time.time()
  tri_dict = nx.triangles(G)
  t1 = time.time()
  print("Time: ", t1-t0)
  total = 0
  for key in tri_dict:
    total += tri_dict[key]
  total = total / 3
  print("Triangles: ",total)

if __name__ == "__main__":
  main()