# Using converter:
`numactl -i all ./converter -bs 32 -rounds 1 -s -m -enc bytepd-amortized -o /ssd1/graphs/tmp/soc-LJ_sym.bytepda ~/inputs/soc-LiveJournal1_sym.adj`
Converts a symmetric adjacencygraph into a bytepd-amortized encoded graph, where
the compression block size is 32.
