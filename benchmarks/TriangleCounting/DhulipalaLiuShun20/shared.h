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
#define EMPTYWTV numeric_limits<size_t>::max()
#define UPDATET1 1
#define UPDATET2 2
#define UPDATET3 3
#define UPDATET4 4
#define UPDATET5 5

#define OLD_EDGE 0
#define NEW_EDGE 1
#define DEL_EDGE 2
#define NO_EDGE -1


namespace gbbs{
namespace DBTGraph{
    using EdgeT = pair<uintE, uintE>;

    struct VtxUpdate{
        uintE id = -1;
        size_t degree = 0; // number of updates
        size_t insert_degree = 0; // number of insertion updates
        size_t insert_low_degree = 0; // number of insertion updates to low
        size_t offset = -1; // offsets in edges

        VtxUpdate(uintE a, size_t o):id(a), offset(o){}
        VtxUpdate(){}
        inline void setDeg(size_t a){degree = a;}
        inline size_t end(){return offset + degree;}
        inline size_t insOffset(){return offset + insert_degree;}

    };

    struct VtxUpdateInsDeg{
        pbbs::sequence<VtxUpdate> vtxNew;
        VtxUpdateInsDeg(pbbs::sequence<VtxUpdate>& a):vtxNew(a){}

        size_t operator ()(size_t i)const {
            return vtxNew[i].insert_degree;
        }
    };

    struct TriangleCounts{ // cache line?
        pbbs::sequence<size_t> c;
        size_t P;
        //c1, c2, c3, c4, c5, c6; c1, c2, c3, c4, c5, c6;
        //TriangleCounts():c1(0),c2(0),c3(0), c4(0), c5(0), c6(0){
        TriangleCounts(){
            P = num_workers();
            c = pbbs::sequence<size_t>(P * 6, (size_t)0);
        }

        //flag \in {1,2,3}
        inline void increment(int flag, size_t val){
            c[6 * worker_id() + flag - 1 ] += val;
        }

        //flag \in {1,2,3}
        inline void decrement(int flag, size_t val){
            c[6 * worker_id() + 2 + flag] += val;
        }

        inline pbbs::sequence<size_t> report(){
            pbbs::sequence<size_t> result = pbbs::sequence<size_t>::no_init(6);
            par_for(0, 6, [&] (size_t i) {
                for(size_t j = 0; j < P; ++j) {
                    result[i] += c[6 * j + i];
                }
            });
            return result;
        }

        inline void clear(){
            par_for(0, 6 * P, [&] (size_t i) {c[i] = 0;});
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
            default:
                cout << "invalid update flag " << flag << endl;
                exit(1);
            }
        }

        inline void cleanUp(){
            c1  = c1 + c2 + c3 -c4 -c5;
        }
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
