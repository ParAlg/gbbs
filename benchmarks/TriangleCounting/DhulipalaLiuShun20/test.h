template <class Graph>
inline void testHaveEdge(DBTGraph::DyGraph<Graph> G){
  using EdgeT = DBTGraph::EdgeT;
  bool a = G.haveEdge(EdgeT(1,2));
  a = G.haveEdge(EdgeT(2,1));
  a = G.haveEdge(EdgeT(0,1));
  a = G.haveEdge(EdgeT(0,4));
}

template <class Graph>
inline void testPreprocessing(DBTGraph::DyGraph<Graph> G){
  using EdgeT = DBTGraph::EdgeT;
  pbbs::sequence<pair<EdgeT, bool>> updates = pbbs::sequence<pair<EdgeT, bool>>::no_init(8);
  updates[0] = make_pair(EdgeT(1,2), true);
  updates[1] = make_pair(EdgeT(3,2), true);
  updates[2] = make_pair(EdgeT(1,2), true);
  updates[3] = make_pair(EdgeT(1,2), false);
  updates[4] = make_pair(EdgeT(1,2), true);
  updates[5] = make_pair(EdgeT(4,2), false);
  updates[6] = make_pair(EdgeT(0,4), false);
  updates[7] = make_pair(EdgeT(1,4), true);

  pbbs::sequence<pair<EdgeT, bool>> updates_final = Preprocessing(G, updates);
}
