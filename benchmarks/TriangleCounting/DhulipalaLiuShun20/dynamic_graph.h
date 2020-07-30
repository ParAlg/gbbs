#pragma once

#include <tuple>
#include "gbbs/gbbs.h"
#include "pbbslib/monoid.h"
// #include "gbbs/pbbslib/sparse_table.h"
#include "sparse_table.h"
#include "set.h"
#include "gbbs/macros.h"

using namespace std;

#define EMPTYV numeric_limits<uintE>::max()
#define EMPTYKVB make_tuple(EMPTYV, 0)
#define EMPTYKV make_tuple(EMPTYV, (SetT *)NULL)
#define UPDATET1 1
#define UPDATET2 2
#define UPDATET3 3
#define UPDATET4 4
#define UPDATET5 5


namespace gbbs{
namespace DBTGraph{
    using EdgeT = pair<uintE, uintE>;

    struct WTV{
        size_t c1, c2, c3;
        WTV():c1(0),c2(0),c3(0){}
        WTV(size_t cc1, size_t cc2, size_t cc3):c1(cc1),c2(cc2),c3(cc3){}
        WTV(size_t flag, size_t val):c1(flag),c2(val),c3(0){}

        // inline size_t getFlag(){return c1;}
        inline size_t getUpdateVal(){return c2;}
        inline void update(const  std::tuple<EdgeT, WTV>& kv){ //(edge key, flag, val, 0)
            size_t flag  = std::get<1>(kv).c1;
            switch(flag) {
            case UPDATET1:
                pbbslib::write_add(&c1, std::get<1>(kv).c2);
                break;
            case UPDATET2:
                pbbslib::write_add(&c2, std::get<1>(kv).c2);
                break;
            case UPDATET3:
                pbbslib::write_add(&c3, std::get<1>(kv).c2);
                break;
            case UPDATET4:
                pbbslib::write_add(&c2, std::get<1>(kv).c2);
                break;
            case UPDATET5:
                pbbslib::write_add(&c3, std::get<1>(kv).c2);
                break;
            default:
                cout << "invalid update flag " << flag << endl;
                exit(1);
            }
            // if(flag == UPDATET1){
            //     pbbslib::write_add(&c1, std::get<1>(kv).c2);
            // }else if(flag == UPDATET2){
            //     pbbslib::write_add(&c2, std::get<1>(kv).c2);
            // }else if(flag == UPDATET3){
            //     pbbslib::write_add(&c3, std::get<1>(kv).c2);
            // }else if(flag == UPDATET4){
            //     pbbslib::write_add(&c2, std::get<1>(kv).c2);
            // }else if(flag == UPDATET5){
            //     pbbslib::write_add(&c3, std::get<1>(kv).c2);
            // }else{
            //     cout << "invalid update flag " << flag << endl;
            //     exit(1);
            // }
        }
    };


    struct vertexHash { //TODO: check
        uint64_t operator ()(const uintE& v) const {return pbbs::hash64_2(v);}
    };

    struct edgeHash { //TODO: check
        uint64_t operator ()(const pair<uintE, uintE>& v) const{
            return pbbs::hash_combine(pbbs::hash64_2(v.first), pbbs::hash64_2(v.second));}
    }; 


    inline pbbslib::sparse_table<uintE, int, vertexHash> make_vertex_set(size_t m, long space_mult=-1) {
        return pbbslib::sparse_table<uintE, int, vertexHash>(m, EMPTYKVB, vertexHash(), space_mult);
    }
    

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

        // bool is_high(size_t k) const { return k > 2*t1;} // 
        // bool is_low(size_t k) const {return !is_high(k);}

        bool is_high_v(uintE v)const{ return status[v];}//is_high(D[v]);}
        bool is_low_v(uintE v)const{return !is_high_v(v);}

        bool must_high(size_t k) const { return k > t2;}
        bool must_low(size_t k) const {return k < t1;}

        bool change_status(uintE v, size_t k) const { return (is_high_v(v) && must_low(k)) || (is_low_v(v) && must_high(k));}

        bool use_block(uintE d)const{return d <= block_size;}
        bool use_block_v(uintE v)const{return use_block(D[v]);}

        inline void insertTop(tableE *tb, uintE u, size_t size, size_t bottom_load ){
            SetT *tbB = new SetT(size, EMPTYKVB, vertexHash(), bottom_load);
            tb->insert(make_tuple(u, tbB));
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
        inline void insertW(uintE u, uintE v, int flag, size_t val = 1){
            // if(u > v) swap(u,v);
            if(u >= v) return; //TODO: only store a wedge once
            T->insert_f(make_tuple(make_pair(u,v), WTV(flag, val)), updateTF());

        }

        inline uintE getEArray(uintE u, size_t i)const {
            return edges[block_size * u + i].first;
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

        void markEdgeArray(uintE u, pbbs::sequence<uintE> vs, int flag){
            size_t degree = D[u];
            parallel_for(0, vs.size(), [&](size_t i) {
                setEArray(u,vs[i],degree+i,flag);
            });
        }

        void markEdgeTable(uintE u, pbbs::sequence<uintE> vs, int flag, bool rev){
            if(is_low_v(u)){
                parallel_for(0, vs.size(), [&](size_t i) { 
                    uintE v = vs[i];
                    if(is_low_v(v)){
                        insertE(LL, u, v, flag);
                        if(rev) insertE(LL, v, u, flag);
                    }else{
                        insertE(LH, u, v, flag);
                        if(rev) insertE(HL, v, u, flag);
                    }
                 });
                
            }else{
                parallel_for(0, vs.size(), [&](size_t i) { 
                    uintE v = vs[i];
                    if(is_low_v(v)){
                        insertE(LH, u, v, flag);
                        if(rev) insertE(HL, v, u, flag);
                    }else{
                        insertE(HH, u, v, flag);
                        if(rev) insertE(HH, v, u, flag);
                    }
                 });
            }
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
                status[i] = D[i] > 2 * t1;
            });

            pbbs::sequence<uintE> vArray = pbbs::sequence<uintE>::no_init(n);
            par_for(0, n, [&] (size_t i) {vArray[i] = i;});
            pbbs::sequence<uintE> highNodes = pbbs::filter(vArray, [=] (size_t i) {return is_high_v(i);});
            size_t lowNum = n - highNodes.size();
            vArray.clear();

            // TODO: can allocate less, count using block nodes
            LL = new tableE(lowNum, EMPTYKV, vertexHash(), 1.0);
            LH = new tableE(lowNum, EMPTYKV, vertexHash(), 1.0);
            HL = new tableE(n-lowNum, EMPTYKV, vertexHash(), 1.0);
            HH = new tableE(n-lowNum, EMPTYKV, vertexHash(), 1.0);
            T  = new tableW((n-lowNum)*(n-lowNum), make_tuple(make_pair(EMPTYV, EMPTYV), WTV()), edgeHash(), 1.0);

            // insert top level keys
            double bottom_load = 1.2;// TUNE
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
                    tableE *tb1 = HL;
                    tableE *tb2 = HH;
                    if(is_low_v(i)){
                        tb1 = LL;
                        tb2 = LH;
                    }
                    if(lowD[i] == 0){ //only has high neighbors
                        insertTop(tb2, i, degree, bottom_load);
                    }else if(lowD[i] < degree){
                        insertTop(tb1, i, lowD[i], bottom_load);
                        insertTop(tb2, i, degree - lowD[i], bottom_load);
                    }else if(lowD[i] == degree){//only has low neighbors
                        insertTop(tb1, i, degree, bottom_load);
                    }else{
                        cout << "dynamic_graph.h wrong size " << lowD[i] << endl;
                        exit(1);
                    }
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