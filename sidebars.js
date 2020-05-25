module.exports = {
  docs: {
    Introduction : ['introduction', 'install',  'run', 'formats', 'inputs', 'python_bindings'],
    Tutorial : ['tutorial/bfs_tutorial'],
    Benchmarks : [
      'overview',
      'definitions',
      {
        type: 'category',
        label: 'Shortest Path Problems',
        items: ['sssp/breadth_first_search'],
      },
      {
        type: 'category',
        label: 'Connectivity Problems',
        items: [
          'connectivity/biconnectivity',
          'connectivity/connectivity',
          'connectivity/low_diameter_decomposition',
          'connectivity/minimum_spanning_forest',
          'connectivity/strongly_connected_components'
        ],
      },
      {
        type: 'category',
        label: 'Covering Problems',
        items: ['covering/mis'],
      },
      {
        type: 'category',
        label: 'Substructure Problems',
        items: ['substructure/kcore'],
      },
      {
        type: 'category',
        label: 'Eigenvector Problems',
        items: ['eigenvector/pagerank'],
      },
    ],
  },
};
