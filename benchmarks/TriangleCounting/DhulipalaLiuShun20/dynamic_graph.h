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

// #define EMPTYV numeric_limits<uintE>::max()
// #define EMPTYKVB make_tuple(EMPTYV, 0)
// #define EMPTYKV make_tuple(EMPTYV, (SetT *)NULL)
// #define UPDATET1 1
// #define UPDATET2 2
// #define UPDATET3 3
// #define UPDATET4 4
// #define UPDATET5 5


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
        pbbs::sequence<size_t> lowD;
        pbbs::sequence<pair<uintE,int>> edges;
        tableE *LL;
        tableE *HH;
        tableE *LH;
        tableE *HL;
        tableW *T;

        bool is_high_v(uintE v)const{ return status[v];}//is_high(D[v]);}
        bool is_low_v(uintE v)const{return !is_high_v(v);}

        bool is_high(size_t k) const { return k > 2*t1;} // 
        bool is_low(size_t k) const {return !is_high(k);}

        bool must_high(size_t k) const { return k > t2;}
        bool must_low(size_t k) const {return k < t1;}

        bool use_block(size_t d)const{return d <= block_size;}
        bool use_block_v(uintE v)const{return use_block(D[v]);}

        inline void insertTop(tableE *tb, uintE u, size_t size, size_t bottom_load = 1.2 ){
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

        inline uintE getEArray(uintE u, size_t i)const {
            return edges[block_size * u + i].first;
        }

        inline uintE getEArrayVal(uintE u, size_t i)const {
            return edges[block_size * u + i].second;
        }

        inline void setEArray(uintE u, uintE v, size_t i, int val){
            edges[block_size * u + i] = make_pair(v,val);
        }

        inline void setEArray(uintE v, size_t k, int val){
            edges[k] = make_pair(v,val);
        }


    public:

        size_t num_vertices() const { return n; }
        size_t num_edges() const { return m; }

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
            if(use_block(degree1)){
                for(size_t i = 0; i < degree1; ++i){
                    if(getEArray(u,i) == v) return true;
                }
                return false;
            }else{
                return tb->find(u, (SetT *)NULL)->contains(v);
            }
        }
        
        //delete if edge found and flag is true
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
            if(use_block(degree1)){
                for(size_t i = 0; i < degree1; ++i){
                    if(getEArray(u,i) == v){ 
                        if(flag) setEArray(u,v,i, DEL_EDGE);
                        return true;}
                }
                return false;
            }else{
                SetT *bottomTb = tb->find(u, (SetT *)NULL);
                if(bottomTb->contains(v)){
                    if(flag) bottomTb->updateSeq(v,DEL_EDGE);
                    return true;
                }else{
                    return false;
                }
                
            }
        }

        /////////////////////// MARK EDGE INSERTION /////////////////////////////////////////////
        // assume there is enough space in array
        void markEdgeArrayInsertion(DBTGraph::VtxUpdate u, pbbs::sequence<pair<EdgeT,bool>> &edgesInsert){
            size_t offset = D[u.id];
            parallel_for(0, u.insert_degree, [&](size_t i) {
                setEArray(u.id, edgesInsert[i].first.second, offset+i, NEW_EDGE);
            });
        }

        // assume in table
        void markEdgeTablesInsertion(DBTGraph::VtxUpdate u, pbbs::sequence<pair<EdgeT,bool>> &edgesInsert, bool resize){
            tableE *tb1 = LL;tableE *tb2 = LH;
            if(is_high_v(u.id)){tb1 = HL;tb2 = HH;}
            if(resize){
                if(lowD[u.id] > 0){
                    SetT *bottomTb1 = tb1->find(u, NULL);
                    bottomTb1->maybe_resize(u.insert_low_degree, lowD[u.id]);
                }
                if(D[u.id]-lowD[u.id] > 0){
                    SetT *bottomTb2 = tb2->find(u, NULL);
                    bottomTb2->maybe_resize(u.insert_degree - u.insert_low_degree,  D[u.id] - lowD[u.id]);
                }
            }
            parallel_for(0, u.insert_degree, [&](size_t i) { 
                uintE v = edgesInsert[i].first.second;
                if(is_low_v(v)){
                    insertE(tb1, u, v, NEW_EDGE);
                }else{
                    insertE(tb2, u, v, NEW_EDGE);
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
        }

        void markEdgeInsertion(DBTGraph::VtxUpdate i, pbbs::range<pair<EdgeT,bool>> &edgesI){
            uintE u = i.id;
            size_t space = D[u] + i.insert_degree;
            if(use_block(space)){ // copy to array
                markEdgeArrayInsertion(i, edgesI);
            }else{ 
                if(use_block_v(u)){// copy from array to table
                    copyArrayToTable(i);
                }
                markEdgeTablesInsertion(i, edgesI, !use_block_v(u));// insert to table
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
        void updateTable(DBTGraph::VtxUpdate w, pbbs::range<pair<EdgeT,bool>> &edgesID){
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
                        par_for(0, H.size(), [&] (size_t j) {
                            uintE v = get<0>(H->table[j]);
                            if(v != H->empty_key && u.id != v){
                                updateTableT(u.id, v, get<1>(H->table[j]), flag);}});                    
                    }
                }
            }
            }
        }
         // pbbs::sequence<DBTGraph::VtxUpdate> &vtxNew, pbbs::sequence<size_t> &vtxMap


        /////////////////////// Count Triangles /////////////////////////////////////////////
        inline void CountTriangles(DBTGraph::VtxUpdate u, DBTGraph::VtxUpdate v, bool flag, TriangleCounts tc){
            

        }


        DyGraph(int t_block_size, Graph& G):block_size(t_block_size){
            n = G.num_vertices();
            m = G.num_edges();
            M  = m + 1; // m already doubled
            t1 = sqrt(M) / 2;
            t2 = 3 * t1;

            D = pbbs::sequence<size_t>(n, [&](size_t i) { return G.get_vertex(i).getOutDegree(); });
            edges = pbbs::sequence<pair<uintE,int>>((size_t)(block_size*n), make_pair(EMPTYV,0));
            lowD = pbbs::sequence<size_t>::no_init(n);
            status = pbbs::sequence<bool>::no_init(n);

            
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
            T  = new tableW((n-lowNum)*(n-lowNum), make_tuple(make_pair(EMPTYV, EMPTYV), WTV()), edgeHash(), 1.0);

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
                    if(is_high_v(u.id)){tb1 = HL;tb2 = HH;}
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