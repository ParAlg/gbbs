            // auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
            //     if(is_low(v_data[v].degree)){flag[]}
            //     // int mark = 1;
            //     // if(is_high(v_data[v].degree)){
            //     //     mark = 2;
            //     // }
            //     // flag[u] = flag[u]|mark; // if low, flag is 1, if high is 2, if both is 3
            // };
            // mapEdges(map_f, true);



// // struct hash_pair {
// //   inline size_t operator () (const std::tuple<uintE, uintE>& t) {
// //     size_t l = std::min(std::get<0>(t), std::get<1>(t));
// //     size_t r = std::max(std::get<0>(t), std::get<1>(t));
// //     size_t key = (l << 32) + r;
// //     return pbbslib::hash64_2(key);
// //   }
// // };

    // struct WTV{
    //     size_t c1, c2, c3;
    //     WTV():c1(0),c2(0),c3(0){}
    //     WTV(size_t cc1, size_t cc2, size_t cc3):c1(cc1),c2(cc2),c3(cc3){}
    //     WTV(size_t flag, size_t val):c1(flag),c2(val),c3(0){}

    //     // inline size_t getFlag(){return c1;}
    //     inline size_t getUpdateVal(){return c2;}
    //     inline void update(const  std::tuple<EdgeT, WTV>& kv){ //(edge key, flag, val, 0)
    //         size_t flag  = std::get<1>(kv).c1;
    //         switch(flag) {
    //         case UPDATET1:
    //             pbbslib::write_add(&c1, std::get<1>(kv).c2);
    //             break;
    //         case UPDATET2:
    //             pbbslib::write_add(&c2, std::get<1>(kv).c2);
    //             break;
    //         case UPDATET3:
    //             pbbslib::write_add(&c3, std::get<1>(kv).c2);
    //             break;
    //         case UPDATET4:
    //             pbbslib::write_add(&c2, std::get<1>(kv).c2);
    //             break;
    //         case UPDATET5:
    //             pbbslib::write_add(&c3, std::get<1>(kv).c2);
    //             break;
    //         default:
    //             cout << "invalid update flag " << flag << endl;
    //             exit(1);
    //         }
    //         // if(flag == UPDATET1){
    //         //     pbbslib::write_add(&c1, std::get<1>(kv).c2);
    //         // }else if(flag == UPDATET2){
    //         //     pbbslib::write_add(&c2, std::get<1>(kv).c2);
    //         // }else if(flag == UPDATET3){
    //         //     pbbslib::write_add(&c3, std::get<1>(kv).c2);
    //         // }else if(flag == UPDATET4){
    //         //     pbbslib::write_add(&c2, std::get<1>(kv).c2);
    //         // }else if(flag == UPDATET5){
    //         //     pbbslib::write_add(&c3, std::get<1>(kv).c2);
    //         // }else{
    //         //     cout << "invalid update flag " << flag << endl;
    //         //     exit(1);
    //         // }
    //     }
    // };



    // struct EdgeT{
    //     uintE first;
    //     uintE second;

    //     EdgeT(uintE a, uintE b):first(a), second(b){}
    //     EdgeT(VtxUpdate a, VtxUpdate b):first(a.id), second(b.id){}
    //     EdgeT(){}
    //     bool operator ==(const EdgeT &b) const{
    //         return b.first ==first && b.second == second;
    //     }

    //     bool operator !=(const EdgeT &b) const{
    //         return b.first !=first || b.second != second;
    //     }

    //     bool operator >(const EdgeT &b) const{
    //         if(first > b.first)return true;
    //         if(first < b.first) return false;
    //         return second > b.second;
    //     }
    // };

  // splits to insertions and deletions, deletions first
  // pbbs::sequence<bool> flag = pbbs::sequence<bool>::no_init(m);

  // par_for(0, m, [&] (size_t i) {
  //     flag[i] = updates_final[i].second;
  // });
  // pair<UpdatesT, size_t> tmp = pbbs::split_two(edges, flag);
  // edges.clear();
  // edges = tmp.first;  
  // m_del = tmp.second;
  // m_ins = m - tmp.second;

template <class EdgeT, class F>
auto sortVertexOutNgh(pbbs::range<pair<EdgeT,bool>> &In, pbbs::range<EdgeT> &Out, F f, flags fl = pbbs::no_flag)
    -> size_t {
  typedef T = pair<EdgeT,bool>;
  size_t n = In.size();
  size_t l = num_blocks(n, pbbs::_block_size);
  sequence<size_t> Sums(l);
  sliced_for(n, pbbs::_block_size,
             [&](size_t i, size_t s, size_t e) {
               size_t c = 0;
               for (size_t j = s; j < e; j++) c += (f(j) == false);
               Sums[i] = c;
             },
             fl);
  size_t m = scan_inplace(Sums.slice(), addm<size_t>());
  sliced_for(n, pbbs::_block_size,
             [&](size_t i, size_t s, size_t e) {
               size_t c0 = Sums[i];
               size_t c1 = s + (m - c0);
               for (size_t j = s; j < e; j++) {
                 if (f(j) == false)
                   assign_uninitialized(Out[c0++], In[j]);
                 else
                   assign_uninitialized(Out[c1++], In[j]);
               }
             },
             fl);
  return m;
}



  // //count D, insert D
  // par_for(0, numVtx, [&] (size_t i) {
  //   size_t s = vtxNew[i].offset;
  //   size_t e = 2*m;
  //   if(i < numVtx) e = vtxNew[i+1].offset;
  //   vtxNew[i].setDeg(e - s);
  //   size_t insert_deg = sortVertexOutNgh(edges.slice(s,e), edges_sorted.slice(s,e), 
  //     [&](const size_t& i) {return !edges[i].second; });

  //   vtxNew[i].insert_degree = insert_deg;
  // });



                bool status_v(uintE u){return is_high_v(u);}
        bool status(size_t k){return is_high(k);}
        bool change_status(uintE v, size_t k) const { return (is_high_v(v) && must_low(k+D[v])) || (is_low_v(v) && must_high(k+D[v]));}

        void moveToNew(tableE *tableFrom, tableE *tableNew, bool changeTo, pbbs::sequence<bool> &newStatus){
            par_for(0, tbFrom->size(), [&] (size_t i) {
                if(tbFrom->table[i] != tbFrom->empty){
                    uintE u = get<0>tbFrom->table[i];
                    if(newStatus[u] == changeTo){
                        tbNew->insert(tbFrom->table[i]);
                    }
                }
            });           
        }
        //new_size = sizeof(to) + fromToTo 
        // TODO: = sizeof(TO) + fromToTo - ToToFrom if resizing
        void rebalanceTop(tableE *tb1,  tableE *tb2, size_t new_size, size_t vtxNum, pbbs::sequence<pair<uintE, size_t>> &vtxNew,
            bool changeTo, pbbs::sequence<bool> &newStatus, 
            pbbs::sequence<size_t> &newMarkD, pbbs::sequence<size_t> &newD, pbbs::sequence<size_t> &newLowD){
            // create tbnew
            tableE *tbFrom1 = HL;
            tableE *tbFrom2 = HH;
            tableE *tbTo1 = LL;
            tableE *tbTo2 = LH;

            if(changeTo){
                swap(tbFrom1, tbTo1);
                swap(tbFrom2, tbTo2);
            }
            tableE *tbNew1 = tbTo1;
            tableE *tbNew2 = tbTo2;
            if(new_size > tbFrom1.size()){
                tbNew1 = new tableE(new_size, EMPTYKV, vertexHash(), 1.0); //HL or LL
                tbNew2 = new tableE(new_size, EMPTYKV, vertexHash(), 1.0); //HH ot LH

                // insert subset of tbFrom to tbnew
                // iterate less?
                cilk_spawn moveToNew(tbFrom1, tbNew1, changeTo, newStatus);
                moveToNew(tbFrom2, tbNew2, changeTo, newStatus);
                cilk_sync;
            }

            // insert subset of tbTo to tbnew, delete from tbTo
            par_for(0, vtxNum, [&] (size_t i) {
                uintE u = vtxNew[i].first;
                if(status[u] != newStatus[u] && newStatus[u] == change){ // move from tbTo to tbnew
                    if(!use_block_v(u)){
                        if(lowD[u] > 0){ //has low neighbors
                            tuple<uintE, SetT*> kv = tbTo1->delVal(u);
                            tbNew1->insert(kv);
                        }
                        if(lowD[u] < D[u]){ //has high neighbors
                            tuple<uintE, SetT*> kv = tbTo2->delVal(u);
                            tbNew2->insert(kv);
                        } 
                    }
                }
            });

            // insert subset of tbFrom to tbTo and insert array to tbFrom and tbNew
            par_for(0, vtxNum, [&] (size_t i) {
                uintE u = vtxNew[i].first;
                if(use_block_v(u) && !use_block(new_degree)){//array to tbFrom or tbNew
                    if(newStatus[u] == changeTo){
                        insertTop(tbNew1, tbNew2, u,is_low(new_degree), newLowD[u] + lowD[u], D[u]+newD[u]);
                    }else{
                        insertTop(tbFrom1, tbFrom2, u,is_low(new_degree), newLowD[u] + lowD[u], D[u]+newD[u]);
                    }
                }else if(status[u] != newStatus[u] && newStatus[u] != changeTo){ // move from tbTo to tbFrom
                    if(!use_block_v(u)){
                        if(lowD[u] > 0){ //has low neighbors
                            tuple<uintE, SetT*> kv = tbFrom1->delVal(u);
                            tbTo1->insert(kv);
                        }
                        if(lowD[u] < D[u]){ //has high neighbors
                            tuple<uintE, SetT*> kv = tbFrom2->delVal(u);
                            tbTo2->insert(kv);
                        } 
                    }
                }
            });


            // insert array to table
        }
        
        //delete if edge found and flag is false
        bool haveEdgeDel (EdgeT e, bool flag) {
            if (e.first >= n || e.second >= n){
                return false;
            }
            uintE u = e.first;
            uintE v = e.second;
            size_t degree1 = D[u];
            size_t degree2 = D[v];
            if(degree1 > degree2) {swap(u,v); swap(degree1, degree2);}
            if(degree1 == 0 || degree2 == 0) return false;
            tableE *tb = LL;
            if(is_high_v(u) && is_high_v(v)){
                tb = HH;
            }else if(is_high_v(v)){
                tb = LH;
            }
            if(use_block_v(u)){
                for(size_t i = 0; i < degree1; ++i){
                    if(getEArray(u,i) == v){ 
                        if(!flag) setEArray(u,v,i, DEL_EDGE);
                        return true;}
                }
                return false;
            }else{
                SetT *bottomTb = tb->find(u, (SetT *)NULL);
                if(bottomTb->contains(v)){
                    if(!flag) bottomTb->updateSeq(v,DEL_EDGE);
                    return true;
                }else{
                    return false;
                }
                
            }
        }



// ======================== major rebalancing ==================================


template <class Graph, class EdgeT>
DBTGraph::DyGraph<Graph> majorRebalancing(DBTGraph::DyGraph<Graph>& DG, pbbs::sequence<pair<EdgeT,bool>> &edges, 
                              pbbs::sequence<DBTGraph::VtxUpdate> vtxNew, size_t new_m){
  //mark deletions in tables
  par_for(0, vtxNew.size(), [&] (size_t i) {
    DG.markEdgeDeletion(vtxNew[i],  edges.slice(vtxNew[i].insOffset(), vtxNew[i].end()));
  });
  
  // //init new tables, update degrees, insert from old graph, insert from array, init T
  DBTGraph::DyGraph<Graph> newDG = DBTGraph::DyGraph<Graph>(DG.get_block_size(), DG.num_vertices(), new_m);
  DG.inherit(newDG, vtxNew);
  DG.clearAfterInherit();// remove old graph

  // count triangles
  return newDG;
}




        void inheritArrays(pbbs::sequence<size_t>& t_D, pbbs::sequence<size_t>& t_lowD, pbbs::sequence<pair<uintE,int>>& t_edges,
                        pbbs::sequence<bool>& t_status, pbbs::sequence<bool>& t_blockStatus,
                        pbbs::sequence<DBTGraph::VtxUpdate> &vtxNew){
            D = pbbs::sequence<size_t>(n, [&](size_t i) { return t_D[i]; });
            lowD = pbbs::sequence<size_t>(n, [&](size_t i) { return t_lowD[i]; });
            status = pbbs::sequence<bool>(n, [&](size_t i) { return t_status[i]; });
            blockStatus = pbbs::sequence<bool>(n, [&](size_t i) { return t_blockStatus[i];});
            edges = t_edges;

            par_for(0, vtxNew.size(), [&] (size_t i) { // now share the same array
                DBTGraph::VtxUpdate u = vtxNew[i];
                updateDegrees(u, false);
                updateStatus(u.id);
                if(t_status(u.id) && use_block_v(u)) packEdgeArrayDeletions(u, t_D[u]);
            });

            auto statusFlag = pbbs::delayed_sequence(status.size(), [&](size_t i)->size_t{
                if(status[i]) return 0;
                return 1;
            });
            auto monoid = pbbslib::addm<size_t>();
            lowNum = pbbs::reduce(statusFlag, monoid);
        }

        void updateStatus(uintE u){
            status[u] = is_high(D[i]); 
            blockStatus[u] = use_block(D[i]);
        }

        void initInsertTopLevelKeys(){
            par_for(0, n, [&] (size_t i) {
                if(!use_block_v(i)){                
                    tableE *tb1 = LL;tableE *tb2 = LH;
                    if(is_high_v(i)){tb1 = HL;tb2 = HH;}
                    insertTop(tb1, tb2, i, lowD[i], D[i]);
                }
            });
        }

        void inheritEdges(uintE u, SetT *tb, DBTGraph::DyGraph<Graph>& newDG){
           par_for(0, tb->size(), [&](size_t i){
                uintE v = get<0>(tb->table[i]);
                if(v != tb->empty_key && get<1>(tb->table[i])!=DEL_EDGE){
                    newDG.insertE(u,v);
                }
            });
        }

        void inheritEdges(uintE u, DBTGraph::DyGraph<Graph>& newDG){
            if(use_block_v(u) && newDG.use_block_v(u)){ //array to array. alredy packed in array init
                return;
            }else if(use_block_v(u)){ //array to table
                for(size_t i = 0; i<D[u]; ++i){
                    if(getEArrayVal(u,i)!=OLD_EDGE){
                        newDG.insertE(u,getEArray(u,i));
                    }
                }
            }else{ 
                tableE *tb1 = LL;tableE *tb2 = LH;
                if(is_high_v(u)){tb1 = HL;tb2 = HH;}
                if((!use_block_v(u)) && newDG.use_block_v(u)){ // table to array
                    size_t offset = 0;
                    if(lowD[u] > 0){
                        offset = newDG.pack_neighbors_helper(tb1->find(u, NULL), u, u*block_size, u*block_size + lowD[u]);
                    }
                    if(lowD[u] <D[u]){
                        newDG.pack_neighbors_helper(tb2->find(u, NULL), u, u*block_size + offset, (u+1)*block_size);
                    }
                }else{ //table to table
                    if(lowD[u] > 0){inheritEdges(u,tb1->find(u,NULL),newDG)}
                    if(lowD[u] <D[u]){inheritEdges(u,tb2->find(u,NULL),newDG)}
                }
            }

        }

        void initEdgeInsertion(DBTGraph::VtxUpdate &u, pbbs::range<pair<EdgeT,bool>> &edgesI, size_t oldDeg){
            if(edgesI.size()==0) return;
            if(use_block_v(u.id)){ // copy to array
                markEdgeArrayInsertion(u, edgesI, OLD_EDGE, oldDeg - u.delDeg());
            }else{ 
                markEdgeTablesInsertion(u, edgesI, false, OLD_EDGE);// insert to table
            }            
        }

        DyGraph(int t_block_size, size_t t_n, size_t t_m):block_size(t_block_size), n(t_n), m(t_m){
            initParams();
        }

        void inherit(DBTGraph::DyGraph<Graph>& newDG, pbbs::sequence<DBTGraph::VtxUpdate> &vtxNew, pbbs::sequence<pair<EdgeT,bool>> &edges){
            newDG.inheritArrays(D,lowD,edges,status,blockStatus, vtxNew);
            newDG.initTables();
            newDG.initInsertTopLevelKeys();
            par_for(0, n, [&](uintE u){
                inheritEdges(u, newDG);
            });
            // insert from arrays
            par_for(0, vtxNew.size(), [&] (size_t i) {
                newDG.initEdgeInsertion(vtxNew[i], edges.slice(vtxNew[i].offset, vtxNew[i].insOffset()), D[vtxNew[i].id]);
            });
            // init T
            par_for(0, n, [&] (size_t i) {
                newDG.insertWedges(i, true); 
            });
        }

        void clearAfterInherit(){
            D.clear();
            lowD.clear();
            status.clear();
            blockStatus.clear();
            clearTableE(LH);
            clearTableE(LL);
            clearTableE(HL);
            clearTableE(HH);
            T->del();
        }


        
./Triangle -s ../../../inputs/rMatGraph_J_5_100