#include <unordered_set>
#include <stack>
#include <random>
#include <cmath>
#include <ctime>
#include <cstdlib>

#include "gbbs/gbbs.h"
#include "gbbs/dynamic_graph_io.h"
#include "gbbs/helpers/sparse_table.h"
#include "benchmarks/TriangleCounting/ShunTangwongsan15/Triangle.h"

namespace gbbs {
using edge_type = std::pair<uintE, uintE>;

struct hashfunc {
    inline size_t operator() (const edge_type& e) const {
        return parlay::hash32(e.first + e.second);
    }
};
using edge_table_type = gbbs::sparse_table<edge_type, bool, hashfunc>;
using endpoint_table_type = gbbs::sparse_table<uintE, bool, std::hash<uintE>>;

double mean(const std::vector<double> &v)
{
    double sum = 0;

    for (auto &each: v)
        sum += each;

    return sum / v.size();
}

double sd(const std::vector<double> &v)
{
    double square_sum_of_difference = 0;
    double mean_var = mean(v);
    auto len = v.size();

    double tmp;
    for (auto &each: v) {
        tmp = each - mean_var;
        square_sum_of_difference += tmp * tmp;
    }

    return std::sqrt(square_sum_of_difference / (len - 1));
}

template<class weight_type>
struct TrianglesPerVertex {

    size_t _n;
    double _error;

    edge_table_type edge_set = gbbs::make_sparse_table<edge_type, bool>(1 << 20,
                        std::make_tuple(std::make_pair(std::numeric_limits<uintE>::max(),
                        std::numeric_limits<uintE>::max()), (bool) 0), hashfunc());
    sequence<endpoint_table_type> adj_list;

    size_t approx_triangles = 0;

    TrianglesPerVertex(size_t n, double error):  _n(n), _error(error) {}

    std::pair<size_t, size_t> CountTrianglesPerVertex(
            std::vector<DynamicEdge<weight_type>> new_edges, size_t start, size_t end) {
            /*BatchDynamicEdges<weight_type>& new_edges, size_t start, size_t end, commandLine& P) {

            auto graph = dynamic_edge_list_to_symmetric_graph(
                new_edges, end);
            auto f = [&] (uintE u, uintE v, uintE w) { };
            auto results = Triangle(graph, f, "kcore", P, true, _error, end);
            return results;*/

        auto contains_a_triangle = gbbs::make_sparse_table<uintE, bool>(1 << 20,
                std::make_tuple(std::numeric_limits<uintE>::max(), (bool) 0),
                std::hash<uintE>{});
        auto edges = sequence<std::pair<uintE, uintE>>(2*(end - start));
        //auto reverse_edges = sequence<std::pair<uintE, uintE>>(end - start);

        parallel_for(start, end, [&](size_t jk) {
            auto from = new_edges[jk].from;
            auto to = new_edges[jk].to;
            auto min_id = std::min(from, to);
            auto max_id = std::max(from, to);

            if (min_id != max_id) {
                edge_set.insert(std::make_tuple(std::make_pair(min_id, max_id), true));
                edges[2*(jk - start)] = std::make_pair(min_id, max_id);
                edges[2*(jk - start) + 1] = std::make_pair(max_id, min_id);
                //reverse_edges[jk - start] = std::make_pair(max_id, min_id);
            }
        });

        auto compare_tup = [&](const std::pair<uintE, uintE>& l,
                const std::pair<uintE, uintE>& r) { return l < r; };
        parlay::sort_inplace(parlay::make_slice(edges), compare_tup);
        //parlay::sort_inplace(parlay::make_slice(reverse_edges), compare_tup);

        auto bool_seq = parlay::delayed_seq<bool>(edges.size() + 1, [&] (size_t i) {
                return (i == 0) || (i == edges.size())
                    || (std::get<0>(edges[i-1]) != std::get<0>(edges[i]));
        });

        /*
        auto bool_seq_reverse = parlay::delayed_seq<bool>(reverse_edges.size() + 1,
                [&] (size_t i) { return (i == 0) || (i == reverse_edges.size())
                    || (std::get<0>(reverse_edges[i-1]) != std::get<0>(reverse_edges[i]));
        });*/

        auto starts = parlay::pack_index(bool_seq);
        //auto reverse_starts = parlay::pack_index(bool_seq_reverse);
        auto counts = sequence<uintE>(2*(end-start), 0);
        auto has_triangle = sequence<uintE>(2*(end-start), 0);
        //(starts.size() - 1 + reverse_starts.size() - 1, 0);

        parallel_for(0, starts.size()-1, [&](size_t jn) {
            //auto length = starts[jn+1] - starts[jn];
            //auto edge_counts = sequence<uintE>(length, 0);
            parallel_for(starts[jn], starts[jn+1], [&] (size_t jk) {
                size_t cur_count = 0;
                auto from = edges[jk].first;
                auto to = edges[jk].second;

                for (size_t ak = jk + 1; ak < starts[jn+1]; ak++) {
                    uintE adj_id = edges[ak].second;

                    auto min_id = std::min(to, adj_id);
                    auto max_id = std::max(to, adj_id);
                    if (edge_set.contains(std::make_pair(min_id, max_id)))
                        cur_count++;
                        //edge_counts[jk - starts[jn]]++;
                }

                size_t num_edges_considered = new_edges.size() - end;
                size_t num_edges = std::min(num_edges_considered, (size_t) 1000000);
                cur_count *= ((num_edges_considered) / _n);
                auto bucket = ceil(log(cur_count)/log(_error));
                counts[jk] = (new_edges.size()/((1.0) * (end - start))) * ceil(pow(_error, bucket));
                if (cur_count > 0)
                    has_triangle[jn] = 1;
            });
            //auto cur_count = parlay::scan_inplace(edge_counts);
                //contains_a_triangle.insert(std::make_tuple(from, true));
        });

    /*
        parallel_for(0, reverse_starts.size()-1, [&](size_t jn) {
            auto from = reverse_edges[reverse_starts[jn]].first;
            if (!contains_a_triangle.contains(from)) {
                auto length = reverse_starts[jn+1] - reverse_starts[jn];
                auto edge_counts = sequence<uintE>(length, 0);
                parallel_for(reverse_starts[jn], reverse_starts[jn+1], [&] (size_t jk) {
                    size_t cur_count = 0;
                    auto to = reverse_edges[jk].second;

                    for (size_t ak = jk + 1; ak < reverse_starts[jn+1]; ak++) {
                        uintE adj_id = reverse_edges[ak].second;

                        auto min_id = std::min(to, adj_id);
                        auto max_id = std::max(to, adj_id);
                        if (edge_set.contains(std::make_pair(min_id, max_id)))
                            edge_counts[jk - reverse_starts[jn]]++;
                    }
                });

                auto cur_count = parlay::scan_inplace(edge_counts);
                auto bucket = ceil(log(cur_count)/log(_error));

                size_t num_edges = new_edges.size()  - end;
                counts[starts.size() - 1 + jn] = (_n/((1.0) * (end - start))) * ceil(pow(_error, bucket));
                if (cur_count > 0) {
                    contains_a_triangle.insert(std::make_tuple(from, true));
                }
            }
        });
        */

        auto total_count = parlay::scan_inplace(counts);
        auto num_has_triangle = parlay::scan_inplace(has_triangle);
        std::cout << "total count: " << total_count << std::endl;
        std::cout << "num has triangle: " << num_has_triangle << std::endl;

        return std::make_pair(total_count, num_has_triangle);
    }
};

template <class W>
inline void ApproximateTriangles (BatchDynamicEdges<W>& batch_edge_list,
        size_t multiplier,
        size_t trials,
        commandLine& P,
        size_t num_vertices,
        double approx_error,
        bool randomTrial = true,
        size_t vertex_sample_size = 1,
        bool randomStream = true){
    auto batch = batch_edge_list.edges;
    auto rng = std::default_random_engine {};

    uintE num_edges = batch.size();
    uintE cutoff = multiplier * ceil(log(num_edges));

    std::vector<double> errors;

    for (size_t trial = 0; trial < trials; trial++) {
        if (randomTrial) {
            std::srand(unsigned(std::time(0)));
            std::random_shuffle(batch.begin(), batch.end());
        }
        uintE batch_size = 1;
        uintE S = 1;
        uintE num_triangles = 0;
        double prob_sampled = 1/sqrt(batch.size());
        uintE num_to_sample = 1;
        size_t skip = 0;
        auto whole_graph = dynamic_edge_list_to_symmetric_graph(
                batch_edge_list, batch.size());
        auto f = [&] (uintE u, uintE v, uintE w) { };
        auto real_results = Triangle(whole_graph, f, "kcore", P);
        auto triangle_real_count = real_results.first;
        size_t approx_count = 0;
        double error = 0;
        TrianglesPerVertex count_struct =
            TrianglesPerVertex<W>((size_t)num_vertices,
                (double) approx_error);

        size_t start = 0;
        for (size_t i = 0; i < batch.size(); i++) {
            if (randomStream) {
                if (num_to_sample == 0) {
                    auto end_size = std::min(i, batch.size());
                    auto graph = dynamic_edge_list_to_symmetric_graph(
                        batch_edge_list,
                        end_size);

                    // Run static triangle counting on graph
                    auto results = Triangle(graph, f, "kcore", P);
                    num_triangles = results.first;

                    if (num_triangles > cutoff) {
                        S = i;
                        approx_count = ceil(num_triangles * ((1.0 * num_edges)/S)
                            * ((1.0 * num_edges - 1)/(S-1)) * ((1.0 * num_edges - 2)/(S-2)));
                        break;
                    }

                    batch_size *= 2;
                    num_to_sample = batch_size;
                } else {
                    num_to_sample--;
                }
            } else {
                if (num_to_sample == 0) {
                    std::cout << "cur start, i: " << start << ", " << i << std::endl;
                    auto outs = count_struct.CountTrianglesPerVertex(batch, start, i + 1);
                    std::cout << "cur start, i: " << start << ", " << i << std::endl;
                    //batch_edge_list);
                    //, 0, i, P);
                    std::cout << "num with triangles: " << outs.second << std::endl;
                    auto triangle_estimate = outs.first;
                    auto non_zero_vertices = outs.second;
                    if (non_zero_vertices / 2 >= vertex_sample_size) {
                        approx_count += triangle_estimate;
                        S += i - start;
                        skip = 1000000000 - (i - start);
                        start = i + skip;
                        batch_size = 1;
                        break;
                    } else {
                        batch_size *= 2;
                    }
                    num_to_sample = batch_size;
                } else if (skip > 0) {
                    skip--;
                } else {
                    num_to_sample--;
                }
            }
        }

        approx_count /= 6;
        if (approx_count == 0)
            approx_count = triangle_real_count;

        error = (1.0 * std::max(triangle_real_count,
            approx_count))/std::min(triangle_real_count, approx_count);

        errors.push_back(error);
        std::cout << "### Trial #: " << trial << std::endl;
        std::cout << "### S: " << S << std::endl;
        std::cout << "### m: " << num_edges << std::endl;
        std::cout << "### estimated triangle count: " << num_triangles << std::endl;
        std::cout << "### Exact Triangle Count: " << triangle_real_count << std::endl;
        std::cout << "### Approx Count: " << approx_count << std::endl;
        std::cout << "### Error: " << error << std::endl;
    }

    auto m_val = mean(errors);
    auto sd_val = sd(errors);
    std::cout << "### Mean: " << m_val << std::endl;
    std::cout << "### Standard Div: " << sd_val << std::endl;
}

}  // namespace gbbs
