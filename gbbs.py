from gbbs_lib import *
import gbbs_lib
import os

def loadGraph(graphPath="",undirected=True, compressed=False):
  if (graphPath == ""):
    print("Expect a non-empty path to the graph input file")
    exit(0)
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


def loadSnap(graphPath="",undirected=True):
  if (graphPath == ""):
    print("Expect a non-empty path to the graph input file")
    exit(0)
  graphName = os.path.splitext(graphPath)[0]+'.adj'
  if (undirected):
    return gbbs_lib.loadSymmetricEdgeListAsGraph(graphPath, graphName)
  else:
    return gbbs_lib.loadAsymmetricEdgeListAsGraph(graphPath, graphName)
  print("Unsupported options")
  exit(0)
