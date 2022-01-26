#include <unordered_set>
#include <stack>
#include <random>
#include <cmath>
#include <ctime>
#include <cstdlib>

#include "gbbs/gbbs.h"
#include "gbbs/dynamic_graph_io.h"
#include "benchmarks/TriangleCounting/ShunTangwongsan15/Triangle.h"

namespace gbbs {

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

template <class W>
inline void ApproximateTriangles (BatchDynamicEdges<W>& batch_edge_list,
        size_t multiplier,
        size_t trials,
        commandLine& P,
        size_t triangles,
        double eps) {
    auto batch = batch_edge_list.edges;

    uintE num_edges = batch.size();

    std::vector<double> errors;

    for (size_t trial = 0; trial < trials; trial++) {
        std::srand(unsigned(std::time(0)));
        std::random_shuffle(batch_edge_list.edges.begin(), batch_edge_list.edges.end());
        uintE num_triangles = 0;
        double r = 1/(sqrt(triangles) * eps);
        uintE S = floor(r * num_edges);
        std::cout << "S: " << S << std::endl;
        size_t approx_count = 1;
        auto f = [&] (uintE u, uintE v, uintE w) { };
        auto whole_graph = dynamic_edge_list_to_symmetric_graph(
            batch_edge_list, batch.size());
        size_t triangle_real_count = Triangle(whole_graph, f, "kcore", P);
        if (S < batch.size()) {
            auto base_graph = dynamic_edge_list_to_symmetric_graph(
                batch_edge_list, S);//std::min(S, batch.size()));
            uintE num_base_triangles = Triangle(base_graph, f, "kcore", P);
            num_triangles += num_base_triangles;
            std::cout << "base_count: " << num_triangles << std::endl;
            auto end_index = S;//std::min(S, batch.size());

            for (size_t i = end_index; i < batch.size(); i++) {
                // Run static triangle counting on the graph with the additional
                // edge
                auto edge_graph = dynamic_edge_list_to_symmetric_graph(
                    batch_edge_list, end_index, i);
                uintE num_edge_triangles = Triangle(edge_graph, f, "kcore", P);
                if (num_edge_triangles > num_base_triangles) {
                    num_triangles += (num_edge_triangles - num_base_triangles);
                    //std::cout << "New Num Triangles: " << num_triangles << std::endl;
                }
            }

            approx_count = ceil((1/(r*r)) * num_triangles);
        } else
            approx_count = triangle_real_count;

        auto error = (1.0 * std::max(triangle_real_count, approx_count))/std::min(triangle_real_count, approx_count);
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
