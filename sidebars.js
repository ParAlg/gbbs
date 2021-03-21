module.exports = {
  docs: {
    Introduction : ['introduction', 'install',  'run', 'formats', 'python_bindings'],
    Tutorial : ['tutorial/bfs_tutorial'],
    Research : ['research', 'contributors'],
    Library : ['library/overview'],
    "Benchmark Implementations" : [
      'benchmarks/overview',
      'benchmarks/definitions',
      {
        type: 'category',
        label: 'Shortest Path Problems',
        items: [
          'benchmarks/sssp/breadth_first_search',
          'benchmarks/sssp/integral_weight_sssp',
          'benchmarks/sssp/positive_weight_sssp',
          'benchmarks/sssp/general_weight_sssp',
          'benchmarks/sssp/ss_widest_path',
          'benchmarks/sssp/ss_betweenness_centrality',
          'benchmarks/sssp/spanner',
        ],
      },
      {
        type: 'category',
        label: 'Connectivity Problems',
        items: [
          'benchmarks/connectivity/biconnectivity',
          'benchmarks/connectivity/connectivity',
          'benchmarks/connectivity/low_diameter_decomposition',
          'benchmarks/connectivity/minimum_spanning_forest',
          'benchmarks/connectivity/spanning_forest',
          'benchmarks/connectivity/strongly_connected_components'
        ],
      },
      {
        type: 'category',
        label: 'Covering Problems',
        items: [
          'benchmarks/covering/maximal_independent_set',
          'benchmarks/covering/maximal_matching',
          'benchmarks/covering/graph_coloring',
          'benchmarks/covering/apx_set_cover',
        ],
      },
      {
        type: 'category',
        label: 'Substructure Problems',
        items: [
          'benchmarks/substructure/triangle_counting',
          'benchmarks/substructure/k_clique_counting',
          'benchmarks/substructure/k_core',
          'benchmarks/substructure/k_truss',
          'benchmarks/substructure/apx_densest_subgraph',
          'benchmarks/substructure/low_outdegree_orientation',
        ],
      },
      {
        type: 'category',
        label: 'Eigenvector Problems',
        items: [
          'benchmarks/eigenvector/pagerank',
          'benchmarks/eigenvector/cosimrank',
        ],
      },
    ],
  },
};
