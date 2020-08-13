#pragma once

#include <tuple>
#include "gbbs/gbbs.h"
#include "pbbslib/monoid.h"
#include "sparse_table.h"
#include "set.h"
#include "gbbs/macros.h"

using namespace std;

#define EMPTYV numeric_limits<uintE>::max()
#define EMPTYVMAP numeric_limits<size_t>::max()
#define EMPTYKVB make_tuple(EMPTYV, 0)
#define EMPTYKV make_tuple(EMPTYV, (SetT *)NULL)
#define EMPTYWTV numeric_limits<size_t>::max()
#define UPDATET1 1
#define UPDATET2 2
#define UPDATET3 3
#define UPDATET4 4
#define UPDATET5 5
#define UPDATECLEAR 6
#define UPDATECLEANUP 7

#define OLD_EDGE 0
#define NEW_EDGE 1
#define DEL_EDGE 2
#define NO_EDGE -1


namespace gbbs{
namespace DBTGraph{
    
    using EdgeT = pair<uintE, uintE>;
    using SymGraph = symmetric_graph<symmetric_vertex, pbbs::empty>;
    using StaticEdgeT = tuple<uintE, pbbs::empty>;

    inline uintE getFirst(pbbs::sequence<pair<EdgeT,bool>> &edges, size_t i){
    return edges[i].first.first;
    }

    inline uintE getSecond(pbbs::sequence<pair<EdgeT,bool>> &edges, size_t i){
    return edges[i].first.second;
    }

    inline uintE getFirst(pbbs::range<pair<EdgeT,bool> *> edges, size_t i){
    return edges[i].first.first;
    }

    inline uintE getSecond(pbbs::range<pair<EdgeT,bool> *> edges, size_t i){
    return edges[i].first.second;
    }

    inline uintE getFirst(pair<EdgeT,bool> e){
    return e.first.first;
    }

    inline uintE getSecond(pair<EdgeT,bool> e){
    return e.first.second;
    }

    inline uintE getFirst(pair<uintE, int> e){
    return e.first;
    }

    inline uintE getSecond(pair<uintE, int> e){
    return e.second;
    }

    template <class W>
    inline uintE getFirst(tuple<uintE,W> e){
    return get<0>(e);
    }

    template <class W>
    inline uintE getSecond(tuple<uintE,W> e){
    return get<1>(e);
    }


    struct VtxUpdate{
        uintE id = -1;
        size_t degree = 0; // number of updates
        size_t insert_degree = 0; // number of insertion updates
        size_t insert_low_degree = 0; // number of insertion updates to low
        size_t offset = -1; // offsets in edges
        size_t delete_low_degree = 0;

        VtxUpdate(uintE a, size_t o):id(a), offset(o){}
        VtxUpdate(uintE a):id(a){}
        inline void setDeg(size_t a){degree = a;}
        inline void setInsDeg(size_t a){insert_degree = a;}
        inline size_t end(){return offset + degree;}
        inline size_t insOffset(){return offset + insert_degree;}
        inline size_t delDeg(){return degree-insert_degree;}
        inline size_t newDeg(size_t oldDeg){return oldDeg + 2*insert_degree - degree;}
        inline size_t newLowDeg(size_t oldDeg){return oldDeg + insert_low_degree - delete_low_degree;}
    };

    struct VtxUpdateInsDeg{
        pbbs::sequence<VtxUpdate> vtxNew;
        VtxUpdateInsDeg(pbbs::sequence<VtxUpdate>& a):vtxNew(a){}

        size_t operator ()(size_t i)const {
            return vtxNew[i].insert_degree;
        }
    };

    template <class SetT>
    struct MakeEdge{
        using T = typename SetT::T; using K = typename SetT::KT;
        uintE u;T *table;K empty_key;
        MakeEdge(uintE uu, T*t_table, K _empty):u(uu), table(t_table), empty_key(_empty){}

        EdgeT operator ()(size_t i)const {
            return EdgeT(get<0>(table[i]),u);
        }
    };

    template <class SetT>
    struct MakeEdgeEntry{
        using T = typename SetT::T; using K = typename SetT::KT;
        uintE u;T *table;K empty_key;
        MakeEdgeEntry(uintE uu, T*t_table, K _empty):u(uu), table(t_table), empty_key(_empty){}

        pair<uintE, int> operator ()(size_t i)const {
            return make_pair(get<0>(table[i]),get<1>(table[i]));
        }
    };

    template <class SetT>
    struct MakeEdgeEntryMajor{
        using T = typename SetT::T; using K = typename SetT::KT;
        uintE u;T *table;K empty_key;
        MakeEdgeEntryMajor(uintE uu, T* t_table, K _empty):u(uu), table(t_table), empty_key(_empty){}

        StaticEdgeT operator ()(size_t i)const {
            if (get<1>(table[i]) == DEL_EDGE) return make_pair(empty_key, pbbslib::empty());
            return make_tuple(get<0>(table[i]), pbbslib::empty());
        }
    };

    template <class SetT>
    struct MakeEdgeLtoH{
        using T = typename SetT::T; using K = typename SetT::KT;
        uintE u;T *table;K empty_key;
        MakeEdgeLtoH(uintE uu, T* t_table, K _empty):u(uu), table(t_table), empty_key(_empty){}

        pair<EdgeT, bool> operator ()(size_t i)const {
            return make_pair(EdgeT(get<0>(table[i]),u), true);
        }
    };

    template <class SetT>
    struct MakeEdgeHtoL{
        using T = typename SetT::T; using K = typename SetT::KT;
        uintE u;T *table;K empty_key;
        MakeEdgeHtoL(uintE uu, T* t_table, K _empty):u(uu), table(t_table), empty_key(_empty){}

        pair<EdgeT, bool> operator ()(size_t i)const {
            return make_pair(EdgeT(get<0>(table[i]),u), false);
        }
    };

    struct VtxRbl{
        uintE id = -1;
        size_t LtoH;
        size_t degree;
        size_t offset = -1; // offsets in Ngh

        VtxRbl(uintE a, size_t o):id(a), offset(o){}
        VtxRbl(uintE a):id(a){}
        inline void setDeg(size_t a){degree = a;}
        inline void setInsDeg(size_t a){LtoH = a;}
        inline size_t newLowDeg(size_t oldDeg){return oldDeg + degree - LtoH - LtoH;}
        inline size_t getHtoL(){return degree-LtoH;}
        inline size_t insOffset(){return offset + LtoH;}
        inline size_t end(){return offset + degree;}

    };

    struct TriangleCounts{ // cache line size 64
        pbbs::sequence<size_t> c; //c1, c2, c3, c4, c5, c6 ... padding to 64; c1, c2, c3, c4, c5, c6; ....padding to 64;
        size_t P;
        static constexpr size_t eltsPerCacheLine = 128 /sizeof(size_t);
        
        //TriangleCounts():c1(0),c2(0),c3(0), c4(0), c5(0), c6(0){
        TriangleCounts(){
            P = num_workers();
            c = pbbs::sequence<size_t>::no_init(P * eltsPerCacheLine);
            clear();
        }

        //flag \in {1,2,3}
        inline void increment(int flag, size_t val){
            c[eltsPerCacheLine * worker_id() + flag - 1 ] += val;
        }

        //flag \in {1,2,3}
        inline void decrement(int flag, size_t val){
            c[eltsPerCacheLine * worker_id() + 2 + flag] += val;
        }

        inline pbbs::sequence<size_t> report(){
            pbbs::sequence<size_t> result = pbbs::sequence<size_t>(6, (size_t)0);
            par_for(0, P, [&] (size_t i){
                for(int j=0;j<6;++j){
                    result[j] += c[eltsPerCacheLine * i + j];
                }
            });
            return result;
        }

        inline void clear(){
            // par_for(0, 6 * P, [&] (size_t i) {c[i] = 0;});
            par_for(0, P, [&] (size_t i){
                for(int j=0;j<6;++j){c[eltsPerCacheLine * i + j] = 0;}
            });
        }
    };

    

    struct WTV{
        size_t c1, c2, c3, c4, c5;
        WTV():c1(0),c2(0),c3(0), c4(0), c5(0){}
        WTV(size_t a):c1(a),c2(a),c3(a), c4(a), c5(a){}

        // WTV(size_t cc1, size_t cc2, size_t cc3):c1(cc1),c2(cc2),c3(cc3){}
        WTV(size_t flag, size_t val):c1(flag),c2(val),c3(0), c4(0), c5(0){}

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
                pbbslib::write_add(&c4, std::get<1>(kv).c2);
                break;
            case UPDATET5:
                pbbslib::write_add(&c5, std::get<1>(kv).c2);
                break;
            case UPDATECLEANUP:
                cleanUp();
                break;
            case UPDATECLEAR:
                if(c1 !=0) c1 = 0;
                break;
            default:
                cout << "invalid update flag " << flag << endl;
                exit(1);
            }
        }

        inline size_t cleanUp(){
            c1  = c1 + c2 + c3 -c4 -c5;
            c2 = 0;
            c3 = 0;
            c4 = 0;
            c5 = 0;
            return c1;
        }

        inline bool operator== (const WTV& j){
            return c1 == j.c1 && c2 == j.c2 && c3 == j.c3 && c4 == j.c4 && c5 == j.c5 ;
        }
        // bool operator == (const WTV& i, const WTV& j){
        //     return i.c1 == j.c1 && i.c2 == j.c2 && i.c3 == j.c3 && i.c4 == j.c4 && i.c5 == j.c5 ;
        // }
    };


    struct vertexHash { //TODO: check
        uint64_t operator ()(const uintE& v) const {return pbbs::hash64_2(v);}
        int cmp(uintE v, uintE b) {return (v > b) ? 1 : ((v == b) ? 0 : -1);}
    };

    struct edgeHash { //TODO: check
        uint64_t operator ()(const EdgeT& v) const{
            return pbbs::hash_combine(pbbs::hash64_2(v.first), pbbs::hash64_2(v.second));}
        int cmp(EdgeT v, EdgeT b) {return (v > b) ? 1 : ((v == b) ? 0 : -1);}

    }; 


    inline pbbslib::sparse_table<uintE, int, vertexHash> make_vertex_set(size_t m, long space_mult=-1) {
        return pbbslib::sparse_table<uintE, int, vertexHash>(m, EMPTYKVB, vertexHash(), space_mult);
    }
    
}}
