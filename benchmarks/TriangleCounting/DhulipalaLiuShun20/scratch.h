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

./Triangle -s ../../../inputs/rMatGraph_J_5_100