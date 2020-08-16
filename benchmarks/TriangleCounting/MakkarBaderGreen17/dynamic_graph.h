#pragma once

#include "gbbs/gbbs.h"
#include "gbbs/macros.h"
#include "pbbslib/monoid.h"
#include "shared.h"
#include "sparse_table.h"
#include <tuple>

#define EMPTYV numeric_limits<uintE>::max()
#define EMPTYINDK numeric_limits<size_t>::max()
#define EMPTYKV make_tuple(EMPTYINDK, EMPTYV)

using namespace std;

namespace gbbs{
namespace DBTGraph{

struct vertex_hash {
  uint64_t operator()(const uintE &v) const { return pbbs::hash64_2(v); }
};

    template <class Graph> //symmetric_graph
    class DyGraph{
    public:
        using vertex = typename Graph::vertex;
        using edge_type = typename Graph::edge_type;
        using tableV = pbbslib::sparse_table<size_t, uintE, vertex_hash>;

        void init_data_structures(size_t n = 0, size_t batch_size = 0) {
          int n_alloc_size = 2 * n;
          int batch_alloc_size = 2 * n;

          cur_edges =
              pbbs::sequence<pbbs::sequence<uintE>>::no_init(n_alloc_size);
          radj_cur_edges = pbbs::sequence<size_t>::no_init(n_alloc_size);
          ind_cur_ids = pbbs::sequence<uintE>::no_init(n_alloc_size);

          ins_edges =
              pbbs::sequence<pbbs::sequence<uintE>>::no_init(batch_alloc_size);
          radj_ins_edges = pbbs::sequence<size_t>::no_init(batch_alloc_size);
          ind_ins_ids = pbbs::sequence<uintE>::no_init(batch_alloc_size);

          del_edges =
              pbbs::sequence<pbbs::sequence<uintE>>::no_init(batch_alloc_size);
          radj_del_edges = pbbs::sequence<size_t>::no_init(batch_alloc_size);
          ind_del_ids = pbbs::sequence<uintE>::no_init(batch_alloc_size);

          id_cur_edges = new tableV(n_alloc_size, EMPTYKV, vertex_hash());
          id_ins_edges = new tableV(n_alloc_size, EMPTYKV, vertex_hash());
          id_del_edges = new tableV(n_alloc_size, EMPTYKV, vertex_hash());

          degrees = pbbs::sequence<size_t>::no_init(n_alloc_size);
          non_zero_deg = 0;
        }

        size_t get_cur_size() { return cur_edges.size(); }

        size_t get_ins_batch_size() { return ins_edges.size(); }

        size_t get_del_batch_size() { return del_edges.size(); }

        // Resize adjacency list method
        void resize_adj_list(size_t new_size, size_t adj_ind,
                             pbbs::sequence<uintE> updates, bool isInsert) {
          size_t j = updates.size() - 1;
          size_t i = cur_edges[adj_ind].size() - 1;

          pbbs::sequence<uintE> *adj_list = cur_edges[adj_ind];

          pbbs::sequence<uintE> new_adj_list =
              pbbs::sequence<size_t>::no_init(new_size);

          if (isInsert) {
            while (i >= 0 && j >= 0) {
              uintE diff = adj_list[i] - updates[j];

              if (diff > 0) {
                new_adj_list[i + j + 1] = adj_list[i];
                i--;
              } else {
                new_adj_list[i + j + 1] = updates[j];
                j--;
              }
            }

            while (j >= 0) {
              new_adj_list[i + j + 1] = updates[j];
            }

            while (i >= 0) {
              new_adj_list[i + j + 1] = adj_list[i];
            }
          } else {
            while (i >= 0 && j >= 0) {
              uintE diff = adj_list[i] - updates[j];
              if (diff == 0) {
                adj_list[i] = -1;
              }

              if (diff >= 0) {
                i--;
              }

              if (diff <= 0) {
                j--;
              }
            }

            i = 0;
            j = 0;
            while (i < new_adj_list.size()) {
              if (adj_list[i] != -1) {
                new_adj_list[j] = adj_list[i];
                j++;
              }
            }
            cur_edges[adj_ind].clear();
            cur_edges[adj_ind] = new_adj_list;
          }
          degrees[adj_ind] = new_size / 2;
        }

        // Merge adjacency list with batch of inserts
        void merge_adj_list(uintE node_id, range<uintE> updates,
                            bool isInsert) {
          if (ind_cur_edges.contains(node_id)) {
            size_t adj_index = ind_cur_edges.find(node_id);
            pbbs::sequence<uintE> *cur_adj_list = cur_edges[adj_index];

            size_t j = updates.size();

            d = degrees[adj_index] - 1;

            if (isInsert) {
              if (d + j > degrees[adj_index]) {
                resize_adj_list(2 * (d + j), adj_index, updates, true);
              } else {
                while (d >= 0 && j >= 0) {
                  uintE diff = *cur_adj_list[d] - updates[j];
                  if (diff > 0) {
                    *cur_adj_list[d + j + 1] = *cur_adj_list[d];
                    i--;
                  } else {
                    *cur_adj_list[d + j + 1] = updates[j];
                    j--;
                  }
                }

                while (j >= 0) {
                  *cur_adj_list[d + j + 1] = insert_updates[j];
                  j--;
                }
              }
              degrees[adj_index] = d + j;
            } else {
              if (d - j < 0.25 * d) {
                resize_adj_list(2 * (d - j), adj_index, updates, false);
              } else {
                while (d >= 0 && j >= 0) {
                  uintE diff = *cur_adj_list[d] - updates[j];
                  if (diff == 0) {
                    *cur_adj_list[d] = -1;
                  }

                  if (diff >= 0) {
                    i--;
                  }

                  if (diff <= 0) {
                    j--;
                  }
                }

                i = 0;
                j = 0;
                while (i < cur_edges[adj_index].size()) {
                  if (adj_list[i] != -1) {
                    cur_edges[adj_index][j] = *cur_adj_list[i];
                    j++;
                  }
                }
              }
              degrees[adj_index] = d - j;
            }

          } else {
              if (non_zero_degree + 1 > cur_edges.size) {
                    resize_edge_list();
              }

              int j = updates.size();

          }
        }

        // Merge the batch of inserts to the current graph to obtain updated
        // graph
        void update_inserts(pbbs::sequence<uintE> inserts,
                            pbbs::sequence<size_t> list_starts,
                            pbbs::sequence<size_t> update_ids) {
          int num_ids = update_ids.ize();
          par_for(0, num_ids - 1, [&](const size_t i) {
            merge_adj_ins_list(
                update_ides[i],
                inserts.slice(list_starts[i], list_starts[i + 1]));
          });
        }

      private:
        // 3 structures for each of the following graphs: the current graph
        // the update graph consisting of edge insertion, and the update
        // graph consisting of edge deletions.

        // Current edges in the graph
        pbbs::sequence<pbbs::sequence<uintE>> cur_edges;
        // Edges that are updates in the insert batch
        pbbs::sequence<pbbs::sequence<uintE>> ins_edges;
        // Edges that are updates in the delete batch
        pbbs::sequence<pbbs::sequence<uintE>> del_edges;

        // Stores the index of each adjacent list to create the reduced
        // adjacency lists
        pbbs::sequence<size_t> radj_cur_edges;
        pbbs::sequence<size_t> radj_ins_edges;
        pbbs::sequence<size_t> radj_del_edges;

        // Vertex sparse table for mapping a vertex ID with the corresponding
        // adjacency list
        tableV *id_cur_edges;
        tableV *id_ins_edges;
        tableV *id_del_edges;

        // Array for mapping indices to vertex IDs
        pbbs::sequence<uintE> ind_cur_ids;
        pbbs::sequence<uintE> ind_ins_ids;
        pbbs::sequence<uintE> ind_del_ids;

        // Degrees list
        pbbs::sequence<size_t> degrees;

        // Set of vertices with nonzero degree
        size_t non_zero_deg;

    }; // DBTGraph
    }  // namespace DBTGraph
