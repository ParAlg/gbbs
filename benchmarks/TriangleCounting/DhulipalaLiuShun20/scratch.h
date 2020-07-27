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

./Triangle -s ../../../inputs/rMatGraph_J_5_100