#pragma once

#include <tuple>
#include "gbbs/gbbs.h"
#include "pbbslib/monoid.h"
#include "sparse_table.h"
#include "set.h"
#include "gbbs/macros.h"
#include "shared.h"

using namespace std;

namespace gbbs{
namespace DBTGraph{


    template <class Graph> //symmetric_graph
    class DyGraph{
    public:
        using vertex = typename Graph::vertex;
        using weight_type = typename Graph::weight_type;
        using edge_type = typename Graph::edge_type;
        using SetT = pbbslib::sparse_table<uintE, int, vertexHash >;
        using tableE = pbbslib::sparse_table<uintE, SetT*, vertexHash >;
        using tableW = pbbslib::sparse_table<EdgeT, WTV, edgeHash>;
    private:
        size_t n;
        size_t m;
        size_t block_size;
        size_t M;
        double t1, t2;
        double threshold;
        size_t lowNum;
        pbbs::sequence<size_t> D;
        pbbs::sequence<bool> status;//true if high
        pbbs::sequence<bool> blockStatus;//true if using block
        pbbs::sequence<size_t> lowD;
        pbbs::sequence<pair<uintE,int>> edges;
        tableE *LL;
        tableE *HH;
        tableE *LH;
        tableE *HL;
        tableW *T;



        bool is_high(size_t k) const { return k > threshold;}
        bool is_low(size_t k) const {return !is_high(k);}

        bool must_high(size_t k) const { return k > t2;}
        bool must_low(size_t k) const {return k < t1;}

        bool use_block(size_t d)const{return d <= block_size;}

        // ------start----- D and lowD are not updated
        bool low_table_empty(DBTGraph::VtxUpdate &i){
            return lowD[i.id] == 0;}

        bool high_table_empty(DBTGraph::VtxUpdate &i){
            return lowD[i.id] == D[i.id];}

        bool new_low_table_empty(DBTGraph::VtxUpdate &i){return get_new_low_degree(i) == 0;}
        bool new_high_table_empty(DBTGraph::VtxUpdate &i){return get_new_low_degree(i) == get_new_degree(i);}

        bool marked_low_table_empty(DBTGraph::VtxUpdate &i){
            return lowD[i.id] == 0 && i.insert_low_degree == 0;}

        bool marked_high_table_empty(DBTGraph::VtxUpdate &i){
            return lowD[i.id] == D[i.id] && i.insert_low_degree == i.insert_degree;}
        // ------end----- D and lowD are not updated

        inline void insertTop(tableE *tb, uintE u, size_t size, double bottom_load = 1.2 ){
            if(size <= 0) return;
            SetT *tbB = new SetT(size, EMPTYKVB, vertexHash(), bottom_load);
            bool suc = tb->insert(make_tuple(u, tbB));
            if(!suc){
                cout << "insertTop duplicated" << endl;
                tbB->clear();
            }
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

                
        inline void insertE(const uintE& u, const uintE& v) {
            if(use_block_v(u)) return; 
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
        }

        struct updateTF { //TODO: check
            void operator () (WTV* v0, const std::tuple<EdgeT, WTV>& kv)const {
                v0->update(kv);
            }
        };

        //need u <= v
        inline void insertW(uintE u, uintE v, int flag, size_t val = 1){
            // if(u > v) swap(u,v);
            if(u >= v) return; //only store a wedge once
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

        inline size_t num_vertices() const { return n; }
        inline size_t num_edges() const { return m; }
        inline size_t get_block_size() const {return block_size;}
        inline size_t num_vertices_low() const {return lowNum;}
        inline void set_vertices_low(size_t a){lowNum = a;}

        inline bool is_high_v(uintE v)const{ return status[v];}//is_high(D[v]);}
        inline bool is_low_v(uintE v)const{return !is_high_v(v);}
        inline bool use_block_v(uintE v)const{return blockStatus[v];}

        inline bool change_status(const DBTGraph::VtxUpdate &v) const { // given id v and new degree k
            size_t new_d = v.newDeg(D[v.id]);
            return (is_high_v(v.id) && must_low(new_d)) || (is_low_v(v.id) && must_high(new_d));}

        inline size_t get_new_degree(DBTGraph::VtxUpdate&u) const {
            return u.newDeg(D[u.id]);}
        
        inline size_t get_degree(uintE u) const {
            return D[u];}

        template <class VTX>
        inline size_t get_new_low_degree(VTX &u) const {
            return u.newLowDeg(lowD[u.id]);}

        inline bool majorRebalance(size_t ins_d, size_t del_d) const {
            size_t new_d = num_edges() + ins_d - del_d;
            return  new_d  < M/4 || new_d  > M;}

        bool haveEdge (EdgeT e) const {
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

        /////////////////////////////// Graph to Array /////////////////////////////////

        // tb is a lower table.
        // put all entries v \in tb into seq_out as (v,u)
        template<class E, class F>
        size_t pack_neighbors_helper(SetT *tb, uintE u, pbbs::range<E *> seq_out) const {
            auto pred = [&](const E& t) { return getFirst(t) != tb->empty_key; };
            auto table_seq = pbbs::delayed_sequence<E, F>(tb->size(), F(u, tb->table, tb->empty_key));
            return pbbslib::filter_out(table_seq, seq_out, pred);
        }

        size_t pack_neighbors_in(SetT *tb, uintE u, size_t s, size_t e) const {
            return pack_neighbors_helper<pair<uintE,int>, MakeEdgeEntry<SetT>>(tb, u, edges.slice(s,e));
        }

        // put the nghs of u into Ngh[ngh_s, ngh_e]. is_low_now is true is u is low before updates (so in LL and/or LH)
        // assume edges array is already updated and packed
        template<class E, class F>
        void get_neighbors_minor(DBTGraph::VtxUpdate &u, pbbs::sequence<E> &Ngh, size_t ngh_s, size_t ngh_e, bool is_low_now) const {
            size_t new_degree = ngh_e - ngh_s;
            if(use_block_v(u.id)){
                for(size_t i  = 0; i<new_degree; ++i){
                    Ngh[ngh_s + i] = make_pair(EdgeT(u.id, getEArray(u.id, i)), is_low_now );
                }
            }else{
                tableE *tb1 = LL; tableE *tb2 = LH;
                if(!is_low_now){tb1 = HL; tb2 = HH;}
                size_t new_low_degree = get_new_low_degree(u);
                size_t tmp = 0;
                if(new_low_degree >0){
                    tmp = pack_neighbors_helper<E, F>(tb1->find(u.id, NULL), u.id, Ngh.slice(ngh_s, new_low_degree));
                    assert(new_low_degree == tmp);}
                if(new_low_degree < new_degree){
                    pack_neighbors_helper<E, F>(tb2->find(u.id, NULL), u.id, Ngh.slice(new_low_degree, ngh_e));}
            }
        }

        void get_neighbors_major(uintE u, pbbs::sequence<StaticEdgeT> &seq_out, size_t offset) const {
            using F = MakeEdgeEntryMajor<SetT>;
            if(use_block_v(u)){
                size_t k = 0;
                for(size_t i = 0; i<D[u]; ++i){ 
                    if(getEArrayVal(u,i)!=DEL_EDGE){
                        seq_out[offset + k] = make_tuple(getEArray(u,i), pbbs::empty());
                        k++;
                    }
                }
            }else{
                tableE *tb1 = LL; tableE *tb2 = LH;
                if(is_high_v(u)){tb1 = HL; tb2 = HH;}
                size_t tmp = 0;
                if(lowD[u] > 0){ 
                    tmp = pack_neighbors_helper<StaticEdgeT, F>(tb1->find(u, NULL), u, seq_out.slice(offset, offset + D[u]));}
                if(lowD[u] < D[u]){
                    pack_neighbors_helper<StaticEdgeT, F>(tb2->find(u, NULL), u, seq_out.slice(tmp, offset + D[u]));}
            }
        }

        /////////////////////// MARK EDGE INSERTION & DELETION /////////////////////////////////////////////
        // assume there is enough space in array
        //  start from offset and end at u.insert_degree
        void markEdgeArrayInsertion(DBTGraph::VtxUpdate &u, pbbs::range<pair<EdgeT,bool>*> &edgesInsert, int val, size_t offset){
            // size_t offset = D[u.id];
            parallel_for(0, u.insert_degree, [&](size_t i) {
                setEArray(u.id, edgesInsert[i].first.second, offset+i, val);
            });
        }

        void markEdgeArrayDeletion(DBTGraph::VtxUpdate &u, pbbs::range<pair<EdgeT,bool>*> &edgesDeletion){
            parallel_for(0, edgesDeletion.size(), [&](size_t j) {
                for(size_t i = 0; i < D[u.id]; ++i) {
                if(getEArray(u.id, i) == edgesDeletion[j].first.second){
                    setEArrayVal(u.id, i, DEL_EDGE);
                }
                }
            });
        }

        // assume in table
        void markEdgeTablesInsertion(DBTGraph::VtxUpdate &u, pbbs::range<pair<EdgeT,bool>*> &edgesInsert, bool resize, int val){
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
                    insertE(tb1, u.id, v, val);
                }else{
                    insertE(tb2, u.id, v, val);
                }
            });

        }

        void copyArrayToTable(DBTGraph::VtxUpdate &u){
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
        
        void markEdgeInsertion(DBTGraph::VtxUpdate &i, pbbs::range<pair<EdgeT,bool> *> edgesI){
            if(edgesI.size()==0) return;
            uintE u = i.id;
            size_t space = D[u] + i.insert_degree;
            if(use_block(space)){ // copy to array
                markEdgeArrayInsertion(i, edgesI, NEW_EDGE, D[i.id]);
            }else{ 
                bool check_resize = true;
                if(use_block_v(u)){// copy from array to table
                    copyArrayToTable(i);
                    check_resize = false;
                }
                markEdgeTablesInsertion(i, edgesI, check_resize, NEW_EDGE);// insert to table
            }
        }

        // if flag is true, update value to val
        // if flag is false, delete value
        void markEdgeTables(DBTGraph::VtxUpdate &u, pbbs::range<pair<EdgeT,bool> *> &edgesM, bool flag, int val){
            tableE *tb1 = LL;tableE *tb2 = LH;
            if(is_high_v(u.id)){tb1 = HL;tb2 = HH;}
            parallel_for(0, edgesM.size(), [&](size_t i) { 
                uintE v = edgesM[i].first.second;
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

        void markEdgeDeletion(DBTGraph::VtxUpdate &i, pbbs::range<pair<EdgeT,bool> *> edgesD){
            if(edgesD.size()==0) return;
            uintE u = i.id;
            if(use_block_v(u)){ // mark in array
                markEdgeArrayDeletion(i, edgesD);
            }else{ 
                markEdgeTables(i, edgesD, true, DEL_EDGE); // mark as old edges
            }
        }

        /////////////////////// update table insertions /////////////////////////////////////////////
        void updateTableT(uintE u, uintE v, int val, bool flag){ // maybe inserting new entries
            if(flag && val == OLD_EDGE){
                if(u > v) swap(u,v); 
                insertW(u, v, UPDATET2, 1);
            }else if(flag && val == NEW_EDGE){
                if(u > v) return; // only insert a wedge once, otherwise if both wedge edges are new inserts, count twice
                insertW(u, v, UPDATET3, 1);
            }else if(!flag && val ==  OLD_EDGE){
                if(u > v) swap(u,v); 
                insertW(u, v, UPDATET4, 1);
            }else if(!flag && val ==  DEL_EDGE){
                if(u > v) return; // only insert a wedge once
                insertW(u, v, UPDATET5, 1);
            }
            // ignore mixed edges
        }
        void updateTableArray(DBTGraph::VtxUpdate &w, uintE u, bool flag){
            par_for(0, D[w.id] + w.insert_degree, [&] (size_t i) { // bruteforce finding high ngh of w
                uintE v = getEArray(w.id, i);
                if(v!=u && is_high_v(v)){
                    updateTableT(u, v, getEArrayVal(w.id, i), flag);
                }
            });
        }

        //edgesID is the insertions and deletions of w in updates
        void updateTable(DBTGraph::VtxUpdate &w, pbbs::range<pair<EdgeT,bool> *> edgesID){
            if (is_low_v(w.id)){ // w is low 
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

        inline void countTrianglesHelper(int val1, int val2, bool flag, TriangleCounts &tc){
            if(val1 == NO_EDGE || val2 == NO_EDGE) return;
            if(flag){// +1 new inserts
                if(val1 == DEL_EDGE || val2 == DEL_EDGE) return;
                tc.increment(val1 + val2 + 1, 1);
            }else{// +1 new deletions
                if(val1 == NEW_EDGE || val2 == NEW_EDGE) return;
                tc.decrement(val1/2 + val2 /2 + 1, 1);
            }
        }

        //tb = XX->find(u)
        inline void countTrianglesHelper(SetT *tb, uintE u, uintE v, bool flag, TriangleCounts &tc){
            par_for(0, tb->size(), [&] (size_t i) {
                uintE w = get<0>(tb->table[i]);
                if(w != tb->empty_key && w != v){
                    int val1 = get<1>(tb->table[i]);
                    int val2 = getEdgeVal(w,v);
                    countTrianglesHelper(val1, val2, flag, tc);
                }
            });

        }

        inline void countTriangles(DBTGraph::VtxUpdate &u, DBTGraph::VtxUpdate &v, bool flag, TriangleCounts &tc){
            if(D[u.id] > D[v.id])swap(u,v);
            if(use_block_v(u.id)){
                par_for(0, D[u.id] + u.insert_degree, [&] (size_t i) {
                    uintE w = getEArray(u.id, i);
                    if(w != v.id){
                    int val1 = getEArrayVal(u.id, i);
                    int val2 = getEdgeVal(w,v.id);
                    countTrianglesHelper(val1, val2, flag, tc);
                    }
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
                if(low_space < space - 1){ //HHH, if only, that's v
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
                        tc.decrement(1, wedges.c1 - wedges.c4 - wedges.c5);
                        tc.decrement(2, wedges.c4);
                        tc.decrement(3, wedges.c5);
                    }
                }
            }

        }

        /////////////////////////////// CLEANUP TABLES /////////////////////////////////
        void packEdgeArrayDeletions(DBTGraph::VtxUpdate &u, size_t e){
            size_t offset = block_size * u.id;
            size_t k= offset;
            for (size_t i = 0; i < e; i++)
                if(edges[offset+i].second != DEL_EDGE) edges[k++] = edges[offset+i];
        }

        // must be called before cleanUpEdgeTablesDeletion because markEdgeArrayInsertion starts from D[u]
        void cleanUpEdgeInsertion(DBTGraph::VtxUpdate &i, pbbs::range<pair<EdgeT,bool> *> edgesI){
            if(edgesI.size()==0) return;
            uintE u = i.id;
            if(use_block_v(u)){ // mark in array
                markEdgeArrayInsertion(i, edgesI, OLD_EDGE, D[i.id]); //before D is updated
            }else{ 
                markEdgeTables(i, edgesI, true, OLD_EDGE); // mark as old edges
            }
        }
        
        // pack deletions in array
        // removes edges from tables
        // to delete empty tables, call updateDegreesDeleteFromTable() after minor rebalancing
        // not packing tables to array now because degrees can change in minor rebalancing
        void cleanUpEdgeDeletion(DBTGraph::VtxUpdate &i, pbbs::range<pair<EdgeT,bool> *> edgesD){
            if(edgesD.size()==0) return;
            uintE u = i.id;
            if(use_block_v(u)){ // pack array
                packEdgeArrayDeletions(i, D[i.id] + i.insert_degree); // deletions will be packed out
            }else{        
                markEdgeTables(i, edgesD, false, OLD_EDGE);// remove edges from tables
            }
        }

        void cleanUpTableT(uintE u, uintE v, bool flag, bool del_flag){
            if(u > v) swap(u,v); 
            if(del_flag){
                WTV wedge = T->find(EdgeT(u,v), WTV(EMPTYWTV));
                if( wedge.c1 != EMPTYWTV && wedge.c1 == 0){
                    T->deleteVal(EdgeT(u,v));
                }
            }else{
                insertW(u, v, UPDATECLEANUP);
            }
        }

        void cleanUpTableArray(DBTGraph::VtxUpdate &w, uintE u, bool flag, bool del_flag){
            par_for(0, D[w.id] + w.insert_degree, [&] (size_t i) { // bruteforce finding high ngh of w
                uintE v = getEArray(w.id, i);
                if(v!=u && is_high_v(v)){
                    cleanUpTableT(u, v, flag, del_flag); //getEArrayVal(u, i), 
                }
            });
        }

        // cleanUp the wedges, let c1 = c1+c2+c3-c4-c5
        // called before tables are cleaned up
        // if del_flag is true, we delete 0 counts wedges, use edgesD if del_flag is true
        void cleanUpTable(DBTGraph::VtxUpdate &w, pbbs::range<pair<EdgeT,bool> *> edgesID, bool del_flag = false){
            if (is_low_v(w.id)){ // w is low and w has high ngh
            par_for(0, edgesID.size(), [&] (size_t i) { // loop over the udpate batch (w,u)
                uintE uid = edgesID[i].first.second;
                if(is_high_v(uid)){               // proceed only when u is high
                    bool flag = edgesID[i].second; // insertion or deletion
                    // if(!(del_flag && flag)){ // continue if not (we want to delete and this is an inserting edge)
                        if(use_block_v(w.id)){           
                            cleanUpTableArray(w ,uid, flag, del_flag );
                        }else{
                            SetT *H = LH->find(w.id, NULL);
                            par_for(0, H->size(), [&] (size_t j) {
                                uintE v = get<0>(H->table[j]);
                                if(v != H->empty_key && uid != v){
                                cleanUpTableT(uid, v, flag, del_flag);}});  //get<1>(H->table[j]),
                        }
                    // }
                }
            });
            }
        }

        /////////////////////////////// Init Graph /////////////////////////////////
        
        void initParams(){
            M  = 2 * m + 1; 
            t1 = sqrt(M) / 2;
            t2 = 3 * t1;
            threshold = 2*t1;
        }

        void initTables(){
            // important: save space in top table for array nodes
            LL = new tableE(lowNum, EMPTYKV, vertexHash(), 1.0);
            LH = new tableE(lowNum, EMPTYKV, vertexHash(), 1.0);
            HL = new tableE(n-lowNum, EMPTYKV, vertexHash(), 1.0);
            HH = new tableE(n-lowNum, EMPTYKV, vertexHash(), 1.0);
            T  = new tableW((size_t)((M/threshold)*(M/threshold)/2 + 1), make_tuple(EdgeT(EMPTYV, EMPTYV), WTV()), edgeHash(), 1.0);        
        }

        DyGraph():n(0), m(0), block_size(0){
        }

        DyGraph(int t_block_size, Graph& G):block_size(t_block_size){
            n = G.num_vertices();
            m = G.num_edges() / 2 ;// edges already doubled
            initParams();

            D = pbbs::sequence<size_t>(n, [&](size_t i) { return G.get_vertex(i).getOutDegree(); });
            edges = pbbs::sequence<pair<uintE,int>>((size_t)(block_size*n), make_pair(EMPTYV,0));
            lowD = pbbs::sequence<size_t>::no_init(n);
            status = pbbs::sequence<bool>(n, [&](size_t i) { return is_high(D[i]); });
            blockStatus = pbbs::sequence<bool>(n, [&](size_t i) { return use_block(D[i]);});

            //compute low degree
            auto monoid = pbbslib::addm<size_t>();
            auto map_f = [&](uintE u, uintE v, const typename Graph::weight_type& wgh) -> size_t {
                if(is_low_v(v)) return 1;
                return 0;
            };
            par_for(0, n, [&] (size_t i) {
                lowD[i] = G.get_vertex(i).template reduceOutNgh<size_t>(i, map_f, monoid);
            });

            pbbs::sequence<uintE> vArray = pbbs::sequence<uintE>::no_init(n);
            par_for(0, n, [&] (size_t i) {vArray[i] = i;});
            pbbs::sequence<uintE> highNodes = pbbs::filter(vArray, [=] (size_t i) {return is_high_v(i);});
            lowNum = n - highNodes.size();
            vArray.clear();

            // important: save space in top table for array nodes
            initTables();

            // insert top level keys
            par_for(0, n, [&] (size_t i) {
                size_t degree = D[i];
                if(use_block(degree)){
                    size_t k = block_size*i; //v_data[i].offset;

                    auto map_f = [&] (const uintE& u, const uintE& v, const typename Graph::weight_type& wgh, size_t ind) {
                        setEArray(v,k + ind,0);
                    };
                    G.get_vertex(i).mapOutNghWithIndex(i, map_f);
                }else{                   
                    tableE *tb1 = LL;tableE *tb2 = LH;
                    if(is_high_v(i)){tb1 = HL;tb2 = HH;}
                    insertTop(tb1, tb2, i, lowD[i], D[i]);
                }
            });

            // insert bottom level
            auto insert_f = [&](const uintE& u, const uintE& v, const typename Graph::weight_type& wgh) {
                insertE(u,v);
            };
            G.mapEdges(insert_f, true);

            // init T
            par_for(0, highNodes.size(), [&] (size_t i) {
                uintE u = highNodes[i];   
                insertWedges(u, false); 
            });

            // cleanup
            highNodes.clear();
        }

        inline void insertWedges(uintE u, bool check_high = false){
            if(check_high && is_low_v(u)) return;
            if(lowD[u] > 0){//(u != HL->empty_key){
            if(use_block_v(u)){
                par_for(0, D[u], [&] (size_t j) {
                    uintE w = getEArray(u,j);//edges[u * block_size + j];
                    if(is_low_v(w)){
                        insertT(u, w);
                    }
                });
            }else{
                    SetT* L = HL->find(u, (SetT*) NULL);
                    par_for(0, L->size(), [&] (size_t j) {
                    uintE w = get<0>(L->table[j]);
                    if(w != L->empty_key && lowD[w] < D[w]-1){
                        insertT(u, w);
                    }
                });                   
            }

            }
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
            // tb->del();
            delete tb;
        }

        void del(){
            D.clear();
            lowD.clear();
            edges.clear();
            status.clear();
            blockStatus.clear();
            clearTableE(LH);
            clearTableE(LL);
            clearTableE(HL);
            clearTableE(HH);
            delete T;
        }

        ~DyGraph(){
            del();
        }

        ///////////////// Minor Rebalance /////////////

        size_t minorRblResizeTop(size_t numHtoL, size_t numLtoH){
            if(num_vertices_low() + numHtoL > LH->size()){
                LH->maybe_resize(numHtoL, num_vertices_low());
                LL->maybe_resize(numHtoL, num_vertices_low());
            }
            size_t numH = num_vertices() - num_vertices_low();
            if(numH + numLtoH > HH->size()){
                HH->maybe_resize(numLtoH, numH);
                HL->maybe_resize(numLtoH, numH);
            }
            return lowNum + numHtoL - numLtoH;
        }

        // prereq: u changes from L to H or H to L, tables are already resized
        // called after moving lower level
        // must insert first, then delete
        // change the status of u when is_delete is true
        void minorRblMoveTopTable(DBTGraph::VtxUpdate &u, bool is_low_now, bool is_delete){
            if(use_block_v(u.id)) return; // not in either tables
            tableE *tb1 = HL;tableE *tb3 = LL; // move from 1 to 3
            tableE *tb2 = HH;tableE *tb4 = LH; // move from 2 to 4
            if(is_low_now){    // currently low
                swap(tb1,tb3);
                swap(tb2,tb4);
            }  
            if(!new_low_table_empty(u)){ 
                if(is_delete){
                    tb1->deleteVal(u.id);
                }else{
                    SetT *L = tb1->find(u.id, NULL);
                    tb3->insert(make_tuple(u.id, L));
                }
            }   
            if(!new_high_table_empty(u)){
                if(is_delete){
                    tb2->deleteVal(u.id);
                }else{                     
                    SetT *H = tb2->find(u.id, NULL);
                    tb4->insert(make_tuple(u.id, H));
                }         
            }
            if(is_delete) status[u.id] = is_low_now;//status true if high after
        }

        // is_low_now is true if v has delta nghbors moving from L to H, otherwise from H to L
        // require: v is not using array, edges is updated and packed, before top tables moves
        void minorRblResizeBottomTable(DBTGraph::VtxUpdate &v, size_t delta, bool is_low_now){
            // if(use_block_v(v.id)) return;
            tableE *tb;
            size_t ne;
            if(is_high_v(v.id) && is_low_now){     
                tb = HH; // we are moving delta entries from HL[v] to HH[v] 
                ne =  get_new_degree(v) - get_new_low_degree(v);  //number of elements in HH[v] now
            }else if(is_high_v(v.id)){             
                tb = HL;
                ne =  get_new_low_degree(v); 
            }else if(is_low_v(v.id) && is_low_now){
                tb = LH;
                ne = get_new_degree(v) - get_new_low_degree(v);  
            }else{                              
                tb = LL;
                ne = get_new_low_degree(v);}

            if(ne == 0){
            insertTop(tb, v.id, delta);
            }else{
            tb->find(v.id, NULL)->maybe_resize(delta, ne);
            }
        }

        // is_low_now is true if v has nghbors Ngh moving from L to H, otherwise from H to L
        // require: edges is updated and packed, tables are already resized, before top tables moves
        void minorRblMoveBottomTable(uintE v, pbbs::range<pair<EdgeT,bool> *> Ngh, bool is_low_now, bool is_delete){
            if(use_block_v(v)) return;
            tableE *tb1 = LL;
            tableE *tb2 = LH;
            if(is_high_v(v)){    
                tb1 = HL; 
                tb2 = HH; // we are moving delta entries from HL[v] to HH[v] 
            }
            if(!is_low_now){swap(tb1,tb2);}   //swap,moving to xL from xH

            SetT *fromSet = tb1->find(v, NULL); // both shoud not be NULL
            SetT *toSet = tb2->find(v, NULL);

            par_for(0, Ngh.size(), [&](const size_t i){
                uintE u = getSecond(Ngh,i);
                if(is_delete){fromSet->deleteVal(u);}
                else{toSet->insert(make_tuple(u,OLD_EDGE));}
            });
        }

        //////////////////////////////////////////// UPDATE T TABLE /////////////////////
        void minorRbldeleteW(uintE u, uintE v){
            if(u == v) return;
            if(u > v) swap(u,v); 
            T->deleteVal(EdgeT(u,v));
        }

        // called before tables and status are are rebalanced, but already updated, u is now H, w is now L (before update)
        // edges are updated and packed based on VtxUpdate 
        // delete all wedges (u,w,v) where v is now H
        // if u.v both change from H to L, update twice
        inline void minorRbldeleteWedgeHelper(uintE u, uintE w, pbbs::sequence<VtxUpdate> &vtxNew,  pbbs::sequence<size_t> &vtxMap){
            if(use_block_v(w)){
                VtxUpdate wobj = VtxUpdate(w);
                if(vtxMap[w] != EMPTYVMAP) wobj = vtxNew[vtxMap[w]];
                par_for(0, get_new_degree(wobj), [&] (size_t k) {
                    uintE v = getEArray(w,k);//edges[w * block_size + k];
                    if(is_high_v(v) && u != v){
                        minorRbldeleteW(u, v);
                        // if(u > v) swap(u,v);
                        // insertW(u, v, UPDATECLEAR);
                    }
                });
            }else{
                SetT* H = LH->find(w,(SetT*) NULL );
                if(H == NULL) return; // we can check if H is NULL beforehand, but makes code messy
                par_for(0, H->size(), [&] (size_t k) {
                    uintE v = get<0>(H->table[k]);
                    if(v != H->empty_key && u != v){
                        minorRbldeleteW(u, v);
                        // if(u > v) swap(u,v);
                        // insertW(u, v, UPDATECLEAR);
                    }
                });
            }

        }

        //called before tables and status are rebalanced, but already updated, u is now H (before update)
        // edges are updated and packed based on VtxUpdate 
        // delete all wedges (u,w,v) where v is now H and w is now L
        // "delete" means change wedge num to 0
        // if u.v both change from H to L, update twice
        void minorRblDeleteWedge(DBTGraph::VtxUpdate &u, pbbs::sequence<VtxUpdate> &vtxNew,  pbbs::sequence<size_t> &vtxMap){
            if(get_new_low_degree(u) > 0){
            if(use_block_v(u.id)){
                par_for(0, get_new_degree(u), [&] (size_t j) {
                    uintE w = getEArray(u.id,j);//edges[u * block_size + j];
                    if(is_low_v(w)){
                        minorRbldeleteWedgeHelper(u.id, w, vtxNew, vtxMap);
                    }
                });
            }else{
                    SetT* L = HL->find(u.id, (SetT*) NULL);
                    par_for(0, L->size(), [&] (size_t j) {
                    uintE w = get<0>(L->table[j]);
                    if(w != L->empty_key){ // we can check if H is NULL here, but makes code messy
                        minorRbldeleteWedgeHelper(u.id, w, vtxNew, vtxMap);
                    }
                });                   
            }

            }
        }


        // u is now H (after update), called after tables AND DEGREES are rebalanced
        void minorRblInsertWedge(DBTGraph::VtxUpdate &u){
            insertWedges(u.id, false);
        }

        //////////////////////////////////////////// Minor Rebalancing CLEANUP /////////////////////

        // update degrees and low degrees from inserts/deletes
        void updateDegrees(DBTGraph::VtxUpdate &u){
            D[u.id] = get_new_degree(u);
            lowD[u.id] = get_new_low_degree(u);
        }

        // update low degrees from rebalancing
        void updateDegrees(DBTGraph::VtxRbl &v){
            lowD[v.id] = get_new_low_degree(v);
        }


        // D and lowD updated to after updates and rebalanced
        // pack table to arrays if new degree is low enough
        // all deleted edges are removed
        // status is not changed, change status in downSizeTablesDeletes
        // free allocated lower table if degree is 0
        void downSizeTables(DBTGraph::VtxUpdate &i){
            uintE u = i.id;
            tableE *tb1 = LL;tableE *tb2 = LH;
            if(is_high_v(u)){tb1 = HL;tb2 = HH;}
            if(use_block(D[u])&&!use_block_v(u)){ // with new degree should block, but is not using it
                size_t offset = u*block_size;
                if(lowD[u] > 0){
                    pack_neighbors_in(tb1->find(u, NULL), u, offset,          offset + lowD[u]);
                }
                if(lowD[u] < D[u]){
                    pack_neighbors_in(tb2->find(u, NULL), u, offset + lowD[u], D[u]);
                }  
            }else if(!use_block(D[u])){ // downsize if w/ new degree should use table
                if(lowD[u]!=0 ){
                    tb1->find(u, NULL)->maybe_resize(lowD[u]);
                }else{
                    SetT *L = tb1->find(u, NULL);
                    if(L != NULL){delete L;} //L->del();free(L);
                }
                if(lowD[u]!=D[u]){
                    tb2->find(u, NULL)->maybe_resize(D[u] - lowD[u]);
                }else{
                    SetT *H = tb2->find(u, NULL);
                    if(H != NULL){delete H;}
                }
            }
        }
        
        // called after degrees are updated&rebalanced and downSizeTables
        // D is now new degrees
        // remove tables that 1) packed to array 2) has zero entry
        void downSizeTablesDeletes(DBTGraph::VtxUpdate &u){
            tableE *tb1 = LL; tableE *tb2 = LH;
            if(is_high_v(u.id)){tb1 = HL; tb2 = HH;}            
            if(use_block(D[u.id])&&!use_block_v(u.id)){ // with new degree should block, but is not using it
                if(lowD[u.id]>0 ){
                    tb1->deleteVal(u.id);
                }
                if(lowD[u.id]<D[u.id]){
                    tb2->deleteVal(u.id);
                }  
                blockStatus[u.id] = true;              
            }else if(!use_block_v(u.id)){ // if table exists it would have zero entry
                if(lowD[u.id]==0 ){
                    tb1->deleteVal(u.id);
                }
                if(lowD[u.id]==D[u.id]){
                    tb2->deleteVal(u.id);
                }    

            }
        }

    };

}
}