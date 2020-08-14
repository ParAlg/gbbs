#pragma once

#include <tuple>
#include "gbbs/pbbslib/sparse_table.h"

// a nested hash table where the keys of two tables have the same type
// DOES NOT WORK YET!!

namespace NestHash{
using namespace std;
using namespace pbbslib;

template <class K, class BV, class KeyHash>
class nested_table {
    using BT = std::tuple<K, BV>; // bottom value
    using V = sparse_table<K, BV, KeyHash>; // top value
    using T = std::tuple<K, V>; // top key value pair

    struct insertBottom {
        BT bottom_kv;
        insertBottom(BT _bottom_kv):bottom_kv(_bottom_kv){}
        void operator ()(std::__tuple_element_t<1,T>* v0, std::tuple<K, V>& kv){
            v0.insert(bottom_kv);
        }
    };

    struct insertBottomS {
        BT bottom_kv;
        insertBottomS(BT _bottom_kv):bottom_kv(_bottom_kv){}
        void operator ()(V* v0, T& kv){
            v0.insert_seq(bottom_kv);
        }
    };

    template<class F>
    struct insertBottomF {
        BT bottom_kv;
        F f;
        insertBottomF(BT _bottom_kv, F _f):bottom_kv(_bottom_kv), f(_f){}
        void operator ()(V* v0, T& kv){
            v0.insert_f(bottom_kv, f);
        }
    };

    public:


        // sparse_table<K, BV, KeyHash>  bottom_table;
        sparse_table<K, V, KeyHash> top_table;
        T* table;
        K empty_key;
        V empty_val = sparse_table<K, BV, KeyHash>();

        size_t size(){
            return top_table.size();
        }

        // static void clearA(T* A, long n, T kv) {
        //     top_table.clearA(A, n, kv);
        // }

        inline size_t hashToRange(size_t h) { return top_table.hashToRange(h); }
        inline size_t firstIndex(K& k) { return top_table.firstIndex(k); }
        inline size_t incrementIndex(size_t h) { return top_table.incrementIndex(h); }

        void del() {
           top_table.del();
        }

        void resize_no_copy(size_t incoming){
            top_table.resize_no_copy(incoming);
        }


        nested_table() {
            top_table = sparse_table<K,V,KeyHash>();
            empty_key = top_table.empty_key;
            table = top_table.table;
        }

        // Size is the maximum number of values the hash table will hold.
        // Overfilling the table could put it into an infinite loop.
        nested_table(size_t _m, BT _empty, KeyHash _key_hash, long inp_space_mult=-1)
        {
            T top_empty = std::make_tuple( std::get<0>(_empty), sparse_table<K,BV,KeyHash>(_m, _empty, _key_hash));
            top_table = sparse_table<K,V,KeyHash>(_m, top_empty, _key_hash, inp_space_mult);
            table = top_table.table;
            empty_key = top_table.empty_key;
        }

        // Size is the maximum number of values the hash table will hold.
        // Overfilling the table could put it into an infinite loop.
        nested_table(size_t _m, BT _empty, KeyHash _key_hash, T* _tab, bool clear=true)
        {
            T top_empty = std::make_tuple( std::get<0>(_empty), sparse_table<K,BV,KeyHash>(_m, _empty, _key_hash));
            top_table = sparse_table<K,V,KeyHash>(_m, top_empty, _key_hash, _tab, clear);
            empty_key = top_table.empty_key;
            table = top_table.table;
        }

        // Pre-condition: k must be present in T.
        inline size_t idx(K k) {
           return top_table.idx(k);
        }

        bool insert(std::tuple<K, V> kv){
            return top_table.insert(kv);
        }



        bool insert(K k, K k2, BV v) {
            T dummy = std::make_tuple(k, empty_val);
            auto insert_bf = [&] (V* v0, const T& tup) {
                v0->insert(std::make_tuple(k2,v));
            };
            // insertBottom(std::make_tuple(k2,v))
            return top_table.insert_f(dummy, insert_bf);
        //     size_t h = firstIndex(k);
        //     while (true) {
        //     if (std::get<0>(table[h]) == empty_key) {
        //         cout << 1 << endl;
        //         if (pbbslib::CAS(&std::get<0>(table[h]), empty_key, k)) {
        // //          std::get<1>(table[h]) = std::get<1>(kv);
        //         cout << 3 << endl;
        //         std::get<1>(table[h]).insert(std::make_tuple(k2,v));
        //         cout << 4 << endl;
        //         return true;
        //         }
        //     }
        //     if (std::get<0>(table[h]) == k) {
        //         cout << 2 << endl;
        //         // f(&std::get<1>(table[h]), kv);
        //         std::get<1>(table[h]).insert(std::make_tuple(k2,v));
        //         return false;
        //     }
        //     h = incrementIndex(h);
        //     }
        //     return false;
        }

        template <class F>
        bool insert_f(std::tuple<K, V> kv, const F& f) {
            return top_table.insert_f(kv, f);
        }

        template <class F>
        bool insert_f(K k, K k2, BV v, const F& f) {
            T dummy = std::make_tuple(k, empty_val);
            return top_table.insert_f(dummy, insertBottomF<F>(make_tuple(k2,v), f));
        //     size_t h = firstIndex(k);
        //     while (true) {
        //     if (std::get<0>(table[h]) == empty_key) {
        //         if (pbbslib::CAS(&std::get<0>(table[h]), empty_key, k)) {
        // //          std::get<1>(table[h]) = std::get<1>(kv);
        //         std::get<1>(table[h]).insert<F>(make_tuple(k2,v), f);
        //         return true;
        //         }
        //     }
        //     if (std::get<0>(table[h]) == k) {
        //         f(&std::get<1>(table[h]), kv);
        //         return false;
        //     }
        //     h = incrementIndex(h);
        //     }
        //     return false;
        }

        bool insert_seq(std::tuple<K, V> kv) {
            return top_table.insert_seq(kv);
        }

        bool insert_seq(K k, K k2, BV v) {
            T dummy = std::make_tuple(k, empty_val);
            return top_table.insert_f(dummy, insertBottomS(make_tuple(k2,v)) );
        }        

        bool insert_check(std::tuple<K, V> kv, bool* abort) {
            return top_table.insert_check(kv, abort);
        }


        bool contains(K k) {
            return top_table.contains(k);
        }

        bool contains(K k, K k2) {
            return top_table.find(k, empty_val).contains(k2);
            // if(top_table.contains(k)){
            //     top_table.find(k, (V)NULL).contains(k2);
            // } else {
            //     return false;
            // }
        }

        V find(K k, V default_value) {
            return top_table.find(k, default_value);
        }

        BV find(K k, K k2, BV default_value) {
            return top_table.find(k, empty_val).find(k2, default_value);
        }

        sequence<T> entries() {//  keys? just pointers?
            return top_table.entries();
        }

        void clear() {// memory leak?
            top_table.clear();
        }

};



template <class K, class BV, class KeyHash>
class array_table {
    using BT = std::tuple<K, BV>; // bottom value
    using V = sparse_table<K, BV, KeyHash>; // top value
    using T = std::tuple<K, V>; // top key value pair

    public:

        sparse_table<K, BV, KeyHash> *top_table;// array of tables
        K empty_key;
        V empty_val = sparse_table<K, BV, KeyHash>();
        size_t m;

        size_t size(){
            return m;
        }

        static void clearA(T* A, long n, T kv) {
            cout << "clearA not implemented" << endl;
            exit(1);
        }

        inline size_t hashToRange(size_t h) { return h;}
        inline size_t firstIndex(K& k) { return k;}
        inline size_t incrementIndex(size_t h) { return h;}

        void del() {
           free(top_table);
        }

        void resize_no_copy(size_t incoming){
          cout << "resize not implemented" << endl;
            exit(1);
        }


        array_table() {
        }

        // Size is the maximum number of values the hash table will hold.
        // Overfilling the table could put it into an infinite loop.
        array_table(size_t _m, BT _empty, KeyHash _key_hash, size_t *_deg, long inp_space_mult=-1)
        {
            top_table = pbbslib::new_array_no_init< sparse_table<K,BV,KeyHash>>(_m);
            par_for(0, _m, pbbslib::kSequentialForThreshold, [&] (size_t i) { 
                top_table[i] = sparse_table<K,BV,KeyHash>(_deg[i], _empty, _key_hash, inp_space_mult);
             });
            empty_key =  std::get<0>(_empty);
        }


        // Pre-condition: k must be present in T.
        inline size_t idx(K k) {
           return k;
        }

        bool insert(K k, K k2, BV v) {
            BT kv = std::make_tuple(k2, v);
            return top_table[k].insert(kv);
        }

        template <class F>
        bool insert_f(K k, K k2, BV v, const F& f) {
            BT kv = std::make_tuple(k2, v);
            return top_table[k].insert_f(kv, f);
        }


        bool insert_seq(K k, K k2, BV v) {
            BT kv = std::make_tuple(k2, v);
            return top_table[k].insert_seq(kv);
        }        

        bool insert_check(std::tuple<K, V> kv, bool* abort) {
            cout << "insert_check not implemented" << endl;
            exit(1);
        }


        bool contains(K k) {
            return k < m;
        }

        bool contains(K k, K k2) {
            return top_table[k].contains(k2);
        }

        BV find(K k, K k2, BV default_value) {
            if (k >= m) return default_value;
            return top_table[k].find(k2, default_value);
        }

        // sequence<T> entries() {//  keys? just pointers?
        //     cout << "entries() not implemented loop table" << endl;
        //     exit(1);
        //     return (sequence<T>)NULL;
        // }

        void clear() {// memory leak?
            top_table.clear();
        }

};


}