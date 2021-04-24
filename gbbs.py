from gbbs_lib import *
import gbbs_lib
import os

def loadGraph(graphPath="",undirected=True, compressed=False, binary=False):
  if (graphPath == ""):
    print("Expect a non-empty path to the graph input file")
    exit(0)
  if (undirected):
    if (compressed):
      return gbbs_lib.readCompressedSymmetricUnweightedGraph(graphPath)
    else:
      return gbbs_lib.readSymmetricUnweightedGraph(graphPath, binary)
  else:
    if (compressed):
      return gbbs_lib.readCompressedAsymmetricUnweightedGraph(graphPath)
    else:
      return gbbs_lib.readAsymmetricUnweightedGraph(graphPath, binary)

def loadFloatGraph(graphPath="", undirected=True):
  if (graphPath == ""):
    print("Expect a non-empty path to the graph input file")
    exit(0)
  if (undirected):
    return gbbs_lib.readSymmetricFloatWeightedGraph(graphPath, binary)
  else:
    return gbbs_lib.readAsymmetricFloatWeightedGraph(graphPath, binary)

def loadSnap(graphPath="", undirected=True):
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

def loadFloatSnap(graphPath="", undirected=True):
  if (graphPath == ""):
    print("Expect a non-empty path to the graph input file")
    exit(0)
  if (undirected):
    return gbbs_lib.loadSymmetricFloatEdgeListAsGraph(graphPath)
  else:
    return gbbs_lib.loadAsymmetricFloatEdgeListAsGraph(graphPath)
  print("Unsupported options")
  exit(0)

def loadFromEdgeList(edges, symmetric=True, weighted=False):
  if (symmetric):
    if (not weighted):
      return gbbs_lib.numpyEdgeListToSymmetricUnweightedGraph(edges)
    else:
      if (edges.dtype == "float64" or edges.dtype == "float32"):
        return gbbs_lib.numpyFloatEdgeListToSymmetricWeightedGraph(edges)
      elif (edges.dtype == "uint64"):
        return gbbs_lib.numpyUintEdgeListToSymmetricWeightedGraph(edges)
      else:
        print("unknown edge type")
        print(edges.dtype)
        exit(0)
  else:
    exit(0)

#def test_numpy_array(np_arr):
#  return gbbs_lib.testNumpyArray(np_arr)
