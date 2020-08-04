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

// no new triangles, all low nodes, triangles2.txt
template <class Graph>
inline void test2(DBTGraph::DyGraph<Graph> G){
  using EdgeT = DBTGraph::EdgeT;
  pbbs::sequence<pair<EdgeT, bool>> updates = pbbs::sequence<pair<EdgeT, bool>>::no_init(8);

  updates[0] = make_pair(EdgeT(7,8), true);//add
  updates[1] = make_pair(EdgeT(3,2), true);
  updates[2] = make_pair(EdgeT(1,2), true);
  updates[3] = make_pair(EdgeT(5,8), false);//remove
  updates[4] = make_pair(EdgeT(1,2), true);
  updates[5] = make_pair(EdgeT(4,2), false);//remove
  updates[6] = make_pair(EdgeT(0,4), false);
  updates[7] = make_pair(EdgeT(1,4), true);//add

}

template <class Graph>
inline void test3(DBTGraph::DyGraph<Graph> G){
  using EdgeT = DBTGraph::EdgeT;
  UpdatesT updates = UpdatesT::no_init(3);
  updates[0] = make_pair(EdgeT(1,2), true);//add
  updates[1] = make_pair(EdgeT(2,4), false);
  updates[2] = make_pair(EdgeT(3,4), false);
}

