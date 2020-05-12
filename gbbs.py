from gbbs_lib import *
import gbbs_lib

def loadGraph(graphPath="",undirected=True, compressed=False):
  if (graphPath == ""):
    print("Expect a non-empty path to the graph input file")
    exit(0)
  print("hello")
  if (undirected):
    if (compressed):
      return gbbs_lib.readCompressedSymmetricUnweightedGraph(graphPath)
    else:
      return gbbs_lib.readSymmetricUnweightedGraph(graphPath)
  else:
    if (compressed):
      return gbbs_lib.readCompressedAsymmetricUnweightedGraph(graphPath)
    else:
      return gbbs_lib.readAsymmetricUnweightedGraph(graphPath)
