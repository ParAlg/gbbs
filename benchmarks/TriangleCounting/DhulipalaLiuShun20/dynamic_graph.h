#pragma once

#include <tuple>
#include "gbbs/gbbs.h"
#include "pbbslib/monoid.h"
// #include "gbbs/pbbslib/sparse_table.h"
#include "sparse_table.h"
#include "set.h"
#include "gbbs/macros.h"
#include "shared.h"

using namespace std;

namespace gbbs{
namespace DBTGraph{


    template <class Graph> //symmetric_graph
    class DyGraph{
        using vertex = typename Graph::vertex;
        using weight_type = typename Graph::weight_type;
        using edge_type = typename Graph::edge_type;
        using SetT = pbbslib::sparse_table<uintE, int, vertexHash >;
        using tableE = pbbslib::sparse_table<uintE, SetT*, vertexHash >;
        using tableW = pbbslib::sparse_table<EdgeT, WTV, edgeHash>;

        size_t n;
        size_t m;
        size_t block_size;
        size_t M;
        double t1, t2;
        pbbs::sequence<size_t> D;
        pbbs::sequence<bool> status;//true if high
        pbbs::sequence<bool> blockStatus;//true if high
        pbbs::sequence<size_t> lowD;
        pbbs::sequence<pair<uintE,int>> edges;
        tableE *LL;
        tableE *HH;
        tableE *LH;
        tableE *HL;
        tableW *T;



        bool is_high(size_t k) const { return k > 2*t1;} // 
        bool is_low(size_t k) const {return !is_high(k);}

        bool must_high(size_t k) const { return k > t2;}
        bool must_low(size_t k) const {return k < t1;}

        bool use_block(size_t d)const{return d <= block_size;}
        bool use_block_v(uintE v)const{return blockStatus[v];}

        inline void insertTop(tableE *tb, uintE u, size_t size, double bottom_load = 1.2 ){
            if(size <= 0) return;
            SetT *tbB = new SetT(size, EMPTYKVB, vertexHash(), bottom_load);
            tb->insert(make_tuple(u, tbB));
        }
        //tb1: *L, tb2: *H
        void insertTop(tableE *tb1, tableE *tb2, uintE i, size_t low_degree, size_t degree){
            double bottom_load = 1.2;
            insertTop(tb1, i, low_degree, bottom_load);// if(low_degree > 0){//has low neighbors// }
            insertTop(tb2, i, degree-low_degree, bottom_load);// if(low_degree < degree){ //has high neighbors// }
            
        }

        inline void insertE(tableE *tb, uintE u, uintE v, int val = 0){
            SetT *tbB = tb->find(u, NULL);
            tbB->insert(make_tuple(v, val));
        }

        struct updateTF { //TODO: check
            void operator () (WTV* v0, const std::tuple<EdgeT, WTV>& kv)const {
                v0->update(kv);
            }
        };

        //need u <= v
        inline void insertW(uintE u, uintE v, int flag, size_t val = 1){
            // if(u > v) swap(u,v);
            if(u >= v) return; //TODO: only store a wedge once
            T->insert_f(make_tuple(EdgeT(u,v), WTV(flag, val)), updateTF());

        }

        inline uintE getEArray(uintE u, size_t i)const { //get ith ngh of u
            return edges[block_size * u + i].first;
        }

        inline int getEArrayVal(uintE u, size_t i)const { //get status of ith ngh of u
            return edges[block_size * u + i].second;
        }

        inline void setEArray(uintE u, uintE v, size_t i, int val){
            edges[block_size * u + i] = make_pair(v,val);
        }

        inline void setEArrayVal(uintE u, size_t i, int val){
            edges[block_size * u + i].second = val;
        }

        inline void setEArray(uintE v, size_t k, int val){
            edges[k] = make_pair(v,val);
        }

        // return OLD_EDGE, NEW_EDGE, DEL_EDGE, or NO_EDGE
        int getEdgeVal (uintE u, uintE v) const{
            if (u >= n || v >= n){
                return false;
            }
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
                    if(getEArray(u,i) == v) return getEArrayVal(u,i);
                }
                return NO_EDGE;
            }else{
                return tb->find(u, (SetT *)NULL)->find(v, NO_EDGE);
            }
        }


    public:

        size_t num_vertices() const { return n; }
        size_t num_edges() const { return m; }

        bool is_high_v(uintE v)const{ return status[v];}//is_high(D[v]);}
        bool is_low_v(uintE v)const{return !is_high_v(v);}

        bool majorRebalance(size_t k){
            return (num_edges()+k) < M/4 || (num_edges()+k) > M;
        }

        bool haveEdge (EdgeT e) const{
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
                    if(getEArray(u,i) == v) return true;
                }
                return false;
            }else{
                return tb->find(u, (SetT *)NULL)->contains(v);
            }
        }

        /////////////////////// MARK EDGE INSERTION /////////////////////////////////////////////
        // assume there is enough space in array
        void markEdgeArrayInsertion(DBTGraph::VtxUpdate u, pbbs::range<pair<EdgeT,bool>*> &edgesInsert, int val){
            size_t offset = D[u.id];
            parallel_for(0, u.insert_degree, [&](size_t i) {
                setEArray(u.id, edgesInsert[i].first.second, offset+i, val);
            });
        }

        void markEdgeArrayDeletion(DBTGraph::VtxUpdate u, pbbs::range<pair<EdgeT,bool>*> &edgesDeletion){
            parallel_for(0, D[u.id], [&](size_t i) {
                parallel_for(0, edgesDeletion.size(), [&](size_t j) {
                if(getEArray(u.id, i) == edgesDeletion[j].first.second){
                    setEArrayVal(u.id, i, DEL_EDGE);
                }
                });
            });
        }

        // assume in table
        void markEdgeTablesInsertion(DBTGraph::VtxUpdate u, pbbs::range<pair<EdgeT,bool>*> &edgesInsert, bool resize){
            tableE *tb1 = LL;tableE *tb2 = LH;
            if(is_high_v(u.id)){tb1 = HL;tb2 = HH;}
            if(resize){
                if(lowD[u.id] > 0){
                    SetT *bottomTb1 = tb1->find(u.id, NULL);
                    bottomTb1->maybe_resize(u.insert_low_degree, lowD[u.id]);
                }
                if(D[u.id]-lowD[u.id] > 0){
                    SetT *bottomTb2 = tb2->find(u.id, NULL);
                    bottomTb2->maybe_resize(u.insert_degree - u.insert_low_degree,  D[u.id] - lowD[u.id]);
                }
            }
            parallel_for(0, u.insert_degree, [&](size_t i) { 
                uintE v = edgesInsert[i].first.second;
                if(is_low_v(v)){
                    insertE(tb1, u.id, v, NEW_EDGE);
                }else{
                    insertE(tb2, u.id, v, NEW_EDGE);
                }
            });

        }

        void copyArrayToTable(DBTGraph::VtxUpdate u){
            tableE *tb1 = LL;tableE *tb2 = LH;
            if(is_high_v(u.id)){tb1 = HL;tb2 = HH;}
            size_t low_space = lowD[u.id] + u.insert_low_degree;
            size_t space = D[u.id] + u.insert_degree;
            insertTop(tb1, u.id, low_space); //// if(low_space > 0){// have low edges }
            insertTop(tb2, u.id, space-low_space);//// if(low_space < space){ }
            
            parallel_for(0, D[u.id], [&](size_t i) {
                uintE v = getEArray(u.id, i);
                if(is_high_v(v)){
                    insertE(tb2,u.id,v,OLD_EDGE);
                }else{
                    insertE(tb1,u.id,v,OLD_EDGE);
                }
            });
            blockStatus[u.id]  = false;
        }

        void markEdgeInsertion(DBTGraph::VtxUpdate i, pbbs::range<pair<EdgeT,bool> *> edgesI){
            if(edgesI.size()==0) return;
            uintE u = i.id;
            size_t space = D[u] + i.insert_degree;
            if(use_block(space)){ // copy to array
                markEdgeArrayInsertion(i, edgesI, NEW_EDGE);
            }else{ 
                bool check_resize = true;
                if(use_block_v(u)){// copy from array to table
                    copyArrayToTable(i);
                    check_resize = false;
                }
                markEdgeTablesInsertion(i, edgesI, check_resize);// insert to table
            }
        }

        // if flag is true, update value to val
        // if flag is false, delete value
        void markEdgeTables(DBTGraph::VtxUpdate u, pbbs::range<pair<EdgeT,bool> *> &edges, bool flag, int val){
            tableE *tb1 = LL;tableE *tb2 = LH;
            if(is_high_v(u.id)){tb1 = HL;tb2 = HH;}
            parallel_for(0, edges.size(), [&](size_t i) { 
                uintE v = edges[i].first.second;
                if(is_low_v(v)){
                    SetT* L = tb1->find(u.id, NULL);
                    if(flag){L->updateSeq(v,val);}
                    else{L->deleteVal(v);}
                }else{
                    SetT* H = tb2->find(u.id, NULL);
                    if(flag){H->updateSeq(v,val);}
                    else{H->deleteVal(v);}
                }
            });
        }

        void markEdgeDeletion(DBTGraph::VtxUpdate i, pbbs::range<pair<EdgeT,bool> *> edgesD){
            if(edgesD.size()==0) return;
            uintE u = i.id;
            if(use_block_v(u)){ // mark in array
                markEdgeArrayDeletion(i, edgesD);
            }else{ 
                markEdgeTables(i, edgesD, true, DEL_EDGE); // mark as old edges
            }
        }

        /////////////////////// update table insertions /////////////////////////////////////////////
        void updateTableT(uintE u, uintE v, int val, bool flag){
            if(u > v) swap(u,v);
            if(flag && val == OLD_EDGE){
                insertW(u, v, 2, 1);
            }else if(flag && val == NEW_EDGE){
                insertW(u, v, 3, 1);
            }else if(!flag && val ==  OLD_EDGE){
                insertW(u, v, 4, 1);
            }else if(!flag && val ==  DEL_EDGE){
                insertW(u, v, 5, 1);
            }
            // ignore mixed edges
        }
        void updateTableArray(DBTGraph::VtxUpdate w, uintE u, bool flag){
            par_for(0, D[w.id] + w.insert_degree, [&] (size_t i) { // bruteforce finding high ngh of w
                uintE v = getEArray(w.id, i);
                if(v!=u && is_high_v(v)){
                    updateTableT(u, v, getEArrayVal(u, i), flag);
                }
            });
        }

        //edgesID is the insertions and deletions of w in updates
        void updateTable(DBTGraph::VtxUpdate w, pbbs::range<pair<EdgeT,bool> *> edgesID){
            if (is_low_v(w.id) && 
               (lowD[w.id]+w.insert_low_degree < D[w.id] + w.insert_degree) ){ // w is low and w has high ngh
            par_for(0, w.degree, [&] (size_t i) { // loop over the udpate batch (w,u)
                uintE uid = edgesID[i].first.second;
                if(is_high_v(uid)){               // proceed only when u is high
                    bool flag = edgesID[i].second; // insertion or deletion
                    if(use_block_v(w.id)){           
                        updateTableArray(w ,uid, flag);
                    }else{
                        SetT *H = LH->find(w.id, NULL);
                        par_for(0, H->size(), [&] (size_t j) {
                            uintE v = get<0>(H->table[j]);
                            if(v != H->empty_key && uid != v){
                                updateTableT(uid, v, get<1>(H->table[j]), flag);}});                    
                    }
                }
            });
            }
        }
         // pbbs::sequence<DBTGraph::VtxUpdate> &vtxNew, pbbs::sequence<size_t> &vtxMap


        /////////////////////// Count Triangles /////////////////////////////////////////////

        inline void countTrianglesHelper(int val1, int val2, bool flag, TriangleCounts tc){
            if(val1 == NO_EDGE || val2 == NO_EDGE) return;
            if(flag){// +1 new inserts
                if(val1 == DEL_EDGE || val2 == DEL_EDGE) return;
                tc.increment(val1 + val2 + 1, 1);
            }else{// +1 new deletions
                if(val1 == NEW_EDGE || val2 == NEW_EDGE) return;
                tc.decrement(val1/2 + val2 /2 + 1, 1);
            }
        }

        inline void countTrianglesHelper(SetT *tb, uintE u, uintE v, bool flag, TriangleCounts tc){
            par_for(0, tb->size(), [&] (size_t i) {
                uintE w = get<0>(tb->table[i]);
                if(w != tb->empty_key && u != w){
                    int val1 = get<1>(tb->table[i]);
                    int val2 = getEdgeVal(w,v);
                    countTrianglesHelper(val1, val2, flag, tc);
                }
            });

        }

        inline void countTriangles(DBTGraph::VtxUpdate u, DBTGraph::VtxUpdate v, bool flag, TriangleCounts tc){
            if(D[u.id] > D[v.id])swap(u,v);
            if(use_block_v(u.id)){
                par_for(0, D[u.id] + u.insert_degree, [&] (size_t i) {
                    uintE w = getEArray(u.id, i);
                    int val1 = getEArrayVal(u.id, i);
                    int val2 = getEdgeVal(w,v.id);
                    countTrianglesHelper(val1, val2, flag, tc);
                }); 
                return ;
            }
            size_t low_space = lowD[u.id] + u.insert_low_degree;
            size_t space = D[u.id] + u.insert_degree;
            if(is_low_v(u.id)){ // at least one low vertex    
                if(low_space > 0){ //LLL or LLH
                    SetT *L = LL->find(u.id, NULL);
                    countTrianglesHelper(L,u.id,v.id,flag, tc);
                }
                if(low_space < space){ //LHL or LHH
                    SetT *H = LH->find(u.id, NULL);
                    countTrianglesHelper(H,u.id,v.id,flag, tc);
                }
            }else{ // both are high vertices
                if(low_space < space){ //HHH
                    SetT *H = HH->find(u.id, NULL);
                    countTrianglesHelper(H,u.id,v.id,flag, tc);
                }
                if(u.id > v.id) swap(u,v);
                WTV wedges = T->find(EdgeT(u.id,v.id), WTV(EMPTYWTV));
                if(wedges.c1!=EMPTYWTV){
                    if(flag){
                        tc.increment(1, wedges.c1);
                        tc.increment(2, wedges.c2);
                        tc.increment(3, wedges.c3);
                    }else{
                        tc.decrement(1, wedges.c1);
                        tc.decrement(2, wedges.c4);
                        tc.decrement(3, wedges.c5);
                    }
                }
            }

        }

        /////////////////////////////// CLEANUP TABLES /////////////////////////////////
        void packEdgeArrayDeletions(DBTGraph::VtxUpdate u){
            size_t offset = block_size * u.id;
            size_t k= offset;
            for (size_t i = 0; i < D[u.id] + u.insert_degree; i++)
                if(edges[offset+i].second != DEL_EDGE) edges[k++] = edges[offset+i];
        }

        void cleanUpEdgeInsertion(DBTGraph::VtxUpdate i, pbbs::range<pair<EdgeT,bool> *> edgesI){
            if(edgesI.size()==0) return;
            uintE u = i.id;
            if(use_block_v(u)){ // mark in array
                markEdgeArrayInsertion(i, edgesI, OLD_EDGE);
            }else{ 
                markEdgeTables(i, edgesI, true, OLD_EDGE); // mark as old edges
            }
        }

        void cleanUpEdgeDeletion(DBTGraph::VtxUpdate i, pbbs::range<pair<EdgeT,bool> *> edgesD){
            if(edgesD.size()==0) return;
            uintE u = i.id;
            if(use_block_v(u)){ // copy to array
                packEdgeArrayDeletions(i); // deletions will be packed out
            }else{ 
                markEdgeTables(i, edgesD, false, OLD_EDGE);// remove edges 
            }
        }

        void cleanUpTableT(uintE u, uintE v){
            if(u > v) swap(u,v);
            size_t h = T->idx(EdgeT(u,v));
            get<1>(T->table[h]).cleanUp();
        }
        void cleanUpTableArray(DBTGraph::VtxUpdate w, uintE u){
            par_for(0, D[w.id] + w.insert_degree, [&] (size_t i) { // bruteforce finding high ngh of w
                uintE v = getEArray(w.id, i);
                if(v!=u && is_high_v(v)){
                    cleanUpTableT(u, v);
                }
            });
        }

        void cleanUpTable(DBTGraph::VtxUpdate w, pbbs::range<pair<EdgeT,bool> *> edgesID){
            if (is_low_v(w.id) && 
               (lowD[w.id]+w.insert_low_degree < D[w.id] + w.insert_degree) ){ // w is low and w has high ngh
            par_for(0, w.degree, [&] (size_t i) { // loop over the udpate batch (w,u)
                uintE uid = edgesID[i].first.second;
                if(is_high_v(uid)){               // proceed only when u is high
                    // bool flag = edgesID[i].second; // insertion or deletion
                    if(use_block_v(w.id)){           
                        cleanUpTableArray(w ,uid);
                    }else{
                        SetT *H = LH->find(w.id, NULL);
                        par_for(0, H->size(), [&] (size_t j) {
                            uintE v = get<0>(H->table[j]);
                            if(v != H->empty_key && uid != v){
                                cleanUpTableT(uid, v);}});                    
                    }
                }
            });
            }
        }


        DyGraph(int t_block_size, Graph& G):block_size(t_block_size){
            n = G.num_vertices();
            m = G.num_edges() / 2 ;// edges already doubled
            M  = 2 * m + 1; 
            t1 = sqrt(M) / 2;
            t2 = 3 * t1;

            D = pbbs::sequence<size_t>(n, [&](size_t i) { return G.get_vertex(i).getOutDegree(); });
            edges = pbbs::sequence<pair<uintE,int>>((size_t)(block_size*n), make_pair(EMPTYV,0));
            lowD = pbbs::sequence<size_t>::no_init(n);
            status = pbbs::sequence<bool>::no_init(n);
            blockStatus = pbbs::sequence<bool>(n, [&](size_t i) { return use_block(D[i]);});

            //compute low degree
            auto monoid = pbbslib::addm<size_t>();
            auto map_f = [&](uintE u, uintE v, const typename Graph::weight_type& wgh) -> size_t {
                if(is_low_v(v)) return 1;
                return 0;
            };
            par_for(0, n, [&] (size_t i) {
                lowD[i] = G.get_vertex(i).template reduceOutNgh<size_t>(i, map_f, monoid);
                status[i] = is_high(D[i]);
            });

            pbbs::sequence<uintE> vArray = pbbs::sequence<uintE>::no_init(n);
            par_for(0, n, [&] (size_t i) {vArray[i] = i;});
            pbbs::sequence<uintE> highNodes = pbbs::filter(vArray, [=] (size_t i) {return is_high_v(i);});
            size_t lowNum = n - highNodes.size();
            vArray.clear();

            // important: save space in top table for array nodes
            LL = new tableE(lowNum, EMPTYKV, vertexHash(), 1.0);
            LH = new tableE(lowNum, EMPTYKV, vertexHash(), 1.0);
            HL = new tableE(n-lowNum, EMPTYKV, vertexHash(), 1.0);
            HH = new tableE(n-lowNum, EMPTYKV, vertexHash(), 1.0);
            T  = new tableW((n-lowNum)*(n-lowNum), make_tuple(EdgeT(EMPTYV, EMPTYV), WTV(EMPTYWTV)), edgeHash(), 1.0);

            // insert top level keys
            par_for(0, n, [&] (size_t i) {
                size_t degree = D[i];
                if(use_block(degree)){
                    size_t k = block_size*i; //v_data[i].offset;

                    auto map_f = [&] (const uintE& u, const uintE& v, const typename Graph::weight_type& wgh, size_t ind) {
                        setEArray(v,k + ind,0);
                    };
                    G.get_vertex(i).mapOutNghWithIndex(i, map_f);

                    // if(is_high_v(i) && lowD[i] != 0){
                    //     insertTop(HL, i, 0, bottom_load);
                    // }
                }else{                   
                    tableE *tb1 = LL;tableE *tb2 = LH;
                    if(is_high_v(i)){tb1 = HL;tb2 = HH;}
                    insertTop(tb1, tb2, i, lowD[i], D[i]);
                    // if(lowD[i] == 0){ //only has high neighbors
                    //     insertTop(tb2, i, degree, bottom_load);
                    // }else if(lowD[i] < degree){
                    //     insertTop(tb1, i, lowD[i], bottom_load);
                    //     insertTop(tb2, i, degree - lowD[i], bottom_load);
                    // }else if(lowD[i] == degree){//only has low neighbors
                    //     insertTop(tb1, i, degree, bottom_load);
                    // }else{
                    //     cout << "dynamic_graph.h wrong size " << lowD[i] << endl;
                    //     exit(1);
                    // }
                }
            });

            // insert bottom level
            auto insert_f = [&](const uintE& u, const uintE& v, const typename Graph::weight_type& wgh) {
                size_t degree = D[u];
                if(use_block(degree)) return; // edges already in array
                if(is_low_v(u)){
                    if(is_low_v(v)){
                        insertE(LL, u,v);
                    }else{
                        insertE(LH, u,v);
                    }
                }else{
                    if(is_low_v(v)){
                        insertE(HL, u,v);
                    }else{
                        insertE(HH, u,v);
                    }
                }
            };
            G.mapEdges(insert_f, true);

            // init T
            // par_for(0, HL->size(), [&] (size_t i) {
            par_for(0, highNodes.size(), [&] (size_t i) {
                // uintE u = get<0>(HL->table[i]);  
                uintE u = highNodes[i];                
                if(lowD[u] > 0){//(u != HL->empty_key){
                if(use_block_v(u)){
                    par_for(0, D[u], [&] (size_t j) {
                        uintE w = getEArray(u,j);//edges[u * block_size + j];
                        if(is_low_v(w)){
                            insertT(u, w);
                        }
                    });
                }else{
                    //  SetT* L = get<1>(HL->table[i]);
                     SetT* L = HL->find(u, (SetT*) NULL);
                     par_for(0, L->size(), [&] (size_t j) {
                        uintE w = get<0>(L->table[j]);
                        if(w != L->empty_key && lowD[w] < D[w]-1){
                            insertT(u, w);
                        }
                    });                   
                }

                }
            });

            // cleanup
            highNodes.clear();
        }
        // u is high, w is low
        inline void insertT(uintE u, uintE w){
            if(use_block_v(w)){
                par_for(0, D[w], [&] (size_t k) {
                    uintE v = getEArray(w,k);//edges[w * block_size + k];
                    if(is_high_v(v)){
                        insertW(u, v, UPDATET1, 1);
                    }
                });
            }else{
                SetT* H = LH->find(w,(SetT*) NULL );
                par_for(0, H->size(), [&] (size_t k) {
                    uintE v = get<0>(H->table[k]);
                    if(v != H->empty_key && u != v){
                        insertW(u, v, UPDATET1, 1);
                    }
                });
            }

        }

        //TODO: better way to clear?
        inline void clearTableE(tableE *tb){
            par_for(0, tb->size(), [&] (size_t i) {
                if(tb->table[i] != tb->empty){
                    std::get<1>(tb->table[i])->del();
                    delete std::get<1>(tb->table[i]);
                }
            });
            tb->clear();
            // delete tb;
        }

        ~DyGraph(){
            D.clear();
            lowD.clear();
            edges.clear();
            clearTableE(LH);
            clearTableE(LL);
            clearTableE(HL);
            clearTableE(HH);
            T->del();
            // delete T;

            //todo: clear up
        }
    };

}
}