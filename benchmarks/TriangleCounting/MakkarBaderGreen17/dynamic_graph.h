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

    // TODO: Need to cleanup pointers and references and also make sure to clear
    // memory.
    template <class Graph> //symmetric_graph
    class DyGraph{
    public:
        using vertex = typename Graph::vertex;
        using edge_type = typename Graph::edge_type;
        using tableV = pbbslib::sparse_table<size_t, uintE, vertex_hash>;

        void init_data_structures(size_t n = 0, size_t batch_size = 0) {
          int n_alloc_size = 2 * n;
          int batch_alloc_size = n;

          cur_edges =
              pbbs::sequence<pbbs::sequence<uintE> *>::no_init(n_alloc_size);
          radj_cur_edges = pbbs::sequence<size_t>::no_init(n_alloc_size);
          ind_cur_ids = pbbs::sequence<uintE>::no_init(n_alloc_size);

          ins_edges =
              pbbs::sequence<pbbs::sequence<uintE>>::no_init(batch_alloc_size);
          ind_ins_ids = pbbs::sequence<uintE>::no_init(batch_alloc_size);

          del_edges =
              pbbs::sequence<pbbs::sequence<uintE>>::no_init(batch_alloc_size);
          ind_del_ids = pbbs::sequence<uintE>::no_init(batch_alloc_size);

          id_cur_edges = new tableV(n_alloc_size, EMPTYKV, vertex_hash());
          id_ins_edges = new tableV(n_alloc_size, EMPTYKV, vertex_hash());
          id_del_edges = new tableV(n_alloc_size, EMPTYKV, vertex_hash());

          degrees = pbbs::sequence<size_t>::no_init(n_alloc_size);
          non_zero_deg = 0;

          cur_triangle_count_ = 0;
        }

        DyGraph(): non_zero_deg(0) {
                init_data_structures();
        }

        DyGraph(pbbs::sequence<EdgeT>& initial_graph_edges, pbbs::sequence<uintE>& initialGraph, pbbs::sequence<uintE>& uniqueIds,
                pbbs::sequence<size_t>& offsets): initial_graph_edges_(initial_graph_edges),
        initial_graph_(initialGraph), initial_ids_(uniqueIds),
        initial_offsets_(offsets){
                size_t n = initial_ids_.size();
                auto max_monoid = pbbslib::maxm<size_t>();

                size_t batch_size = pbbs::scan_inplace(offsets.slice(), max_monoid);

                init_data_structures(n, batch_size);
                create_update_graph(initial_graph_, initial_offsets_, initial_ids_, true);
                update_inserts(initial_graph_, initial_offsets_, initial_ids_);

                cur_triangle_count_ = countTriangles(initial_graph_edges_, true);
        }

        size_t get_cur_size() { return non_zero_deg; }
        size_t get_ins_batch_size() { return ins_edges.size(); }
        size_t get_del_batch_size() { return del_edges.size(); }
        size_t get_cur_triangle_count() { return cur_triangle_count_;}

        // Resize adjacency list method
        void resize_adj_list(size_t new_size, size_t adj_ind,
                             pbbs::sequence<uintE> updates, bool isInsert) {
          size_t j = updates.size() - 1;
          size_t i = cur_edges[adj_ind].size() - 1;

          uintE id = ind_cur_ids[adj_ind];

          size_t reduce_index;
          if (isInsert)
              reduce_index = new_size;
            else
              reduce_index = 0;

          pbbs::sequence<uintE> new_adj_list =
              pbbs::sequence<uintE>::no_init(new_size);

          if (isInsert) {
            while (i >= 0 && j >= 0) {
              uintE diff = cur_edges[adj_ind][i] - updates[j];

              if (diff > 0) {
                new_adj_list[i + j + 1] = cur_edges[adj_ind][i];
                if (id > cur_edges[adj_ind][i])
                    reduce_index = i + j + 1;
                i--;
              } else {
                new_adj_list[i + j + 1] = updates[j];
                if (id > updates[j])
                    reduce_index = i + j + 1;
                j--;
              }
            }

            while (j >= 0) {
              new_adj_list[i + j + 1] = updates[j];
                if (id > updates[j])
                    reduce_index = i + j + 1;
            }

            while (i >= 0) {
              new_adj_list[i + j + 1] = cur_edges[adj_ind][i];
                if (id > cur_edges[adj_ind][i])
                    reduce_index = i + j + 1;
            }
          } else {
            while (i >= 0 && j >= 0) {
              uintE diff = cur_edges[adj_ind][i] - updates[j];
              if (diff == 0) {
                cur_edges[adj_ind][i] = NULL;
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
              if (cur_edges[adj_ind][i] != NULL) {
                new_adj_list[j] = cur_edges[adj_ind][i];
                if (id <= cur_edges[adj_ind][i])
                    reduce_index = j;
                j++;
              }
            }
            reduce_index++;
            cur_edges[adj_ind].clear();
            cur_edges[adj_ind] = new_adj_list;
          }
          degrees[adj_ind] = new_size / 2;
        }

        void resize_edge_list() {
          size_t s = cur_edges.size();
          pbbs::sequence<pbbs::sequence<uintE> *> new_edges_list =
              pbbs::sequence<pbbs::sequence<uintE> *>::no_init(2 * s);
          par_for(0, 2 * s,
                  [&](const size_t i) { new_edges_list[i] = &cur_edges[i]; });
          cur_edges.clear();
          cur_edges = new_edges_list;

          // Also need to resize the table for storing node IDs here.
            id_cur_edges->clear();
            tableV *newTable = new tableV(2 * s, EMPTYKV, vertex_hash());

            par_for(0, non_zero_deg, [&] (const size_t i) {
                uintE id = ind_cur_ids[i];
                newTable->insert(std::tuple<uintE, size_t>(id, i));
            });

            id_cur_edges = newTable;
        }

        void create_new_adj_list(uintE node_id, size_t j,
                                 pbbs::sequence<uintE> updates) {
          cur_edges[non_zero_deg] = pbbs::sequence<uintE>::no_init(2 * j);
          par_for(0, j, [&](const size_t i) {
            cur_edges[non_zero_deg][i] = updates[i];
          });
        }

        // Merge adjacency list with batch of inserts
        // TODO: make sure the reduce adj list methods are propagated throughout
        // even when doing the resizings. Currently it is not.
        void merge_adj_list(uintE node_id, pbbs::sequence<uintE> updates,
                            bool isInsert) {
          if (id_cur_edges->contains(node_id)) {
            size_t adj_index = id_cur_edges->find(node_id, 0);
            pbbs::sequence<uintE> cur_adj_list = cur_edges[adj_index];

            size_t j = updates.size();

            size_t d = degrees[adj_index] - 1;

            // Find the first index where the ID in the adjacency list is
            // greater than the node id
            size_t reduce_index;
            if (isInsert)
              reduce_index = d + j;
            else
                reduce_index = 0;

            if (isInsert) {
              if (d + j > degrees[adj_index]) {
                resize_adj_list(2 * (d + j), adj_index, updates, true);
              } else {
                while (d >= 0 && j >= 0) {
                  uintE diff = cur_adj_list[d] - updates[j];
                  if (diff > 0) {
                    cur_adj_list[d + j + 1] = cur_adj_list[d];
                    if (cur_adj_list[d] > node_id)
                      reduce_index = d + j + 1;
                    d--;
                  } else {
                    cur_adj_list[d + j + 1] = updates[j];
                    if (updates[j] > node_id)
                      reduce_index = d + j + 1;
                    j--;
                  }
                }

                while (j >= 0) {
                  cur_adj_list[d + j + 1] = updates[j];
                  if (updates[j] > node_id)
                    reduce_index = d + j + 1;
                  j--;
                }
              }
              degrees[adj_index] = d + j;
            } else {
              if (d - j < 0.25 * d) {
                resize_adj_list(2 * (d - j), adj_index, updates, false);
              } else {
                while (d >= 0 && j >= 0) {
                  uintE diff = cur_adj_list[d] - updates[j];
                  if (diff == 0) {
                    cur_adj_list[d] = NULL;
                  }

                  if (diff >= 0) {
                    d--;
                  }

                  if (diff <= 0) {
                    j--;
                  }
                }

                size_t i = 0;
                j = 0;
                while (i < cur_edges[adj_index].size()) {
                  if (cur_adj_list[i] != NULL) {
                    cur_edges[adj_index][j] = cur_adj_list[i];
                    if (cur_adj_list[i] >= node_id)
                      reduce_index = j;
                    j++;
                  }
                }
                reduce_index++;
              }
              degrees[adj_index] = d - j;
            }

          } else {
            if (non_zero_deg + 1 > cur_edges.size()) {
              resize_edge_list();
            }

            int j = updates.size();
            create_new_adj_list(node_id, j, updates);
            non_zero_deg++;
          }
        }

        // Merge the batch of inserts to the current graph to obtain updated
        // graph
        // list_starts should have size update_ids.size() + 1
        // TODO: check the above
        void update_inserts(pbbs::sequence<uintE> inserts,
                            pbbs::sequence<size_t> list_starts,
                            pbbs::sequence<size_t> update_ids) {
          size_t num_ids = update_ids.size();
          par_for(0, num_ids - 1, [&](const size_t i) {
            merge_adj_list(update_ids[i],
                           inserts.slice(list_starts[i], list_starts[i + 1]),
                           true);
          });
        }

        // Merge the batch of deletes to the current graph to obtain updated
        // graph
        void update_deletes(pbbs::sequence<uintE> deletes,
                            pbbs::sequence<size_t> list_starts,
                            pbbs::sequence<size_t> update_ids) {
          size_t num_ids = update_ids.size();
          par_for(0, num_ids - 1, [&](const size_t i) {
            merge_adj_list(update_ids[i],
                           deletes.slice(list_starts[i], list_starts[i + 1]),
                           true);
          });
        }

        // Create the update graph from insertions and/or deletion batched
        // updates
        // TODO: Make sure offsets is length +1 greater than the lengths of the
        // other lists
        void create_update_graph(pbbs::sequence<uintE> updates,
                                 pbbs::sequence<size_t> offsets,
                                 pbbs::sequence<uintE> ids, bool isInsert) {
          if (isInsert) {
            par_for(0, ids.size(), [&](const size_t i) {
              uintE u = ids[i];
              ins_edges[i] = updates.slice(offsets[i], offsets[i + 1]);
              ind_ins_ids[i] = u;
              id_ins_edges[u] = i;
            });
          } else {
            par_for(0, ids.size(), [&](const size_t i) {
              uintE u = ids[i];
              del_edges[i] = updates.slice(offsets[i], offsets[i + 1]);
              ind_del_ids[i] = u;
              id_del_edges[u] = i;
            });
          }
        }

        // Counts the number of intersections between adjacency lists
        size_t countIntersection(pbbs::sequence<uintE> &adj_1,
                                 pbbs::sequence<uintE> &adj_2, bool use_reduced,
                                 uintE u = 0) {
          size_t count = 0;
          size_t i, j = 0;

          size_t ind_u = ind_cur_ids[u];

          size_t a = adj_1.size();
          size_t b = adj_2.size();

          if (use_reduced) {
            pbbs::sequence<uintE> first = adj_1.slice(0, radj_cur_edges[ind_u]);
          } else {
            pbbs::sequence<uintE> first = adj_1;
          }

          while (i < a && b < j) {
            size_t diff = adj_1[i] - adj_2[j];
            if (diff == 0)
              count += 1;
            else if (diff <= 0)
              i++;
            else if (diff >= 0)
              j++;
          }

          return count;
        }

        size_t countTriangles(pbbs::sequence<EdgeT> updates, bool isInsert) {
          pbbs::sequence<size_t> updates_triangles_count =
              pbbs::sequence<size_t>::no_init(updates.size());
          // First count all the update intersections in the current graph
          par_for(0, updates_triangles_count.size(), [&](const size_t i) {
            uintE id_u = ind_cur_ids[updates[i].first];
            uintE id_v = ind_cur_ids[updates[i].second];

            size_t u = id_cur_edges->find(id_u, 0);
            size_t v = id_cur_edges->find(id_v, 0);

            size_t num_intersects =
                countIntersection(cur_edges[u].slice(0, degrees[u]),
                        cur_edges[v].slice(0, degrees[v]), false);
            updates_triangles_count[i] = num_intersects;
          });

          // Then count all the update intersections between current graph
          // and the update graph
          par_for(0, updates_triangles_count.size(), [&](const size_t i) {
            uintE id_u = ind_cur_ids[updates[i].first];
            uintE id_v = ind_cur_ids[updates[i].second];

            size_t u = id_cur_edges->find(id_u, 0);

            pbbs::sequence<uintE> adj_list;
            if (isInsert) {
              size_t v = id_ins_edges->find(id_v, 0);
              adj_list = ins_edges[v];
            } else {
              size_t v = id_del_edges->find(id_v, 0);
              adj_list = del_edges[v];
            }

            size_t num_intersects =
                countIntersection(cur_edges[u].slice(0, degrees[u]), adj_list, true, u);
            if (isInsert)
              updates_triangles_count[i] -= num_intersects;
            else
              updates_triangles_count[i] += num_intersects;
          });

          // Need a separate sequence because we need to divide this total by
          // three Because we are not truncating adjacency lists
          pbbs::sequence<size_t> update_graph_counts =
              pbbs::sequence<size_t>::no_init(updates.size());

          // Finally count all the intersections in the update graph
          par_for(0, updates_triangles_count.size(), [&](const size_t i) {
            uintE id_u = ind_cur_ids[updates[i].first];
            uintE id_v = ind_cur_ids[updates[i].second];

            pbbs::sequence<uintE> adj_list_u;
            pbbs::sequence<uintE> adj_list_v;
            if (isInsert) {
              size_t u = id_ins_edges->find(id_u, 0);
              size_t v = id_ins_edges->find(id_v, 0);
              adj_list_u = ins_edges[u];
              adj_list_v = ins_edges[v];
            } else {
              size_t v = id_del_edges->find(id_v, 0);
              size_t u = id_del_edges->find(id_u, 0);
              adj_list_u = del_edges[u];
              adj_list_v = del_edges[v];
            }

            size_t num_intersections =
                countIntersection(adj_list_u, adj_list_v, false);
            update_graph_counts[i] = num_intersections;
          });

          // Parallel scan to get the total
          auto monoid = pbbslib::addm<size_t>();
          size_t t_sum =
              pbbs::scan_inplace(updates_triangles_count.slice(), monoid);

          size_t u_sum =
              pbbs::scan_inplace(update_graph_counts.slice(), monoid);
          return t_sum + u_sum / 3;
        }

        // Clean up all unnecessary structures
        void cleanup() {
          ins_edges.clear();
          del_edges.clear();
          id_ins_edges->clear();
          id_del_edges->clear();
          ind_ins_ids.clear();
          ind_del_ids.clear();
        }
        // Clean up initialization structures
        void cleanup_initial() {
          initial_graph_.clear();
          initial_graph_edges_.clear();
          initial_ids_.clear();
          initial_offsets_.clear();
        }

      private:
        // 3 structures for each of the following graphs: the current graph
        // the update graph consisting of edge insertion, and the update
        // graph consisting of edge deletions.

        // Current edges in the graph
        pbbs::sequence<pbbs::sequence<uintE> *> cur_edges;
        // Edges that are updates in the insert batch
        pbbs::sequence<pbbs::sequence<uintE>> ins_edges;
        // Edges that are updates in the delete batch
        pbbs::sequence<pbbs::sequence<uintE>> del_edges;

        // Stores the index of each adjacent list to create the reduced
        // adjacency lists
        pbbs::sequence<size_t> radj_cur_edges;

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

        // Initial graph information
        pbbs::sequence<uintE> initial_graph_;
        pbbs::sequence<EdgeT> initial_graph_edges_;
        pbbs::sequence<uintE> initial_ids_;
        pbbs::sequence<size_t> initial_offsets_;

        // Current triangle count
        size_t cur_triangle_count_;

    }; // class end
    }  // namespace DBTGraph
    }  // namespace gbbs
