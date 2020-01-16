#!/usr/bin/python3.5

import sys
import os

graphs = []
#graphs = ["/ssd1/graphs/bench_experiments/soc-LiveJournal1_sym.adj", "/ssd1/graphs/bench_experiments/com-orkut.ungraph.adj", "/ssd0/graphs/bench_experiments/usa_road.adj", "/ssd1/graphs/bench_experiments/twitter_sym.adj", "/ssd1/graphs/bench_experiments/friendster_sym.adj"]
#c_graphs = ["/ssd0/graphs/bench_experiments/clueweb_sym.bytepda", "/ssd0/graphs/bench_experiments/hyperlink2014_sym.bytepda", "/ssd1/graphs/bench_experiments/hyperlink2012_sym.bytepda"]
c_graphs = ["/ssd1/graphs/bench_experiments/hyperlink2012_sym.bytepda"]
binary_dir = "mains"
output_dir = "data"
rounds = 3

os.makedirs(output_dir, exist_ok=True)


binaries=["bfscc", "gbbscc", "jayanti_bfs", "jayanti_kout", "jayanti_ldd", "jayanti_nosample", "label_propagation", "liutarjan_nosample", "liutarjan_bfs", "liutarjan_ldd", "liutarjan_kout", "shiloach_vishkin", "unite_bfs", "unite_early_bfs", "unite_early_kout", "unite_early_ldd", "unite_early_nosample", "unite_kout", "unite_ldd", "unite_nd_bfs", "unite_nd_kout", "unite_nd_ldd", "unite_nd_nosample", "unite_nosample", "unite_rem_cas_bfs", "unite_rem_cas_kout", "unite_rem_cas_ldd", "unite_rem_cas_nosample", "unite_rem_lock_bfs", "unite_rem_lock_kout", "unite_rem_lock_ldd", "unite_rem_lock_nosample"]


bad_binaries=["bfscc", "jayanti_nosample", "label_propagation", "liutarjan_nosample", "shiloach_vishkin", "unite_early_nosample", "unite_nd_nosample", "unite_nosample"]


def drop_caches():
  os.system('sync; echo 1 | sudo tee /proc/sys/vm/drop_caches')

def run_cmd(binary, rounds, graph, output_file):
  cmd = 'sudo numactl -i all ' + './' + binary_dir + '/' + binary \
            + ' -s -m ' \
            + '-r ' + str(rounds) + ' ' \
            + graph + ' ' \
            + '>> ' + output_file
  print(cmd)
  os.system(cmd)

def run_c_cmd(binary, rounds, graph, output_file):
  cmd = 'sudo numactl -i all ' + './' + binary_dir + '/' + binary \
            + ' -s -c -m ' \
            + '-r ' + str(rounds) + ' ' \
            + graph + ' ' \
            + '>> ' + output_file
  print(cmd)
  os.system(cmd)

for graph in graphs:
  name = os.path.basename(graph)
  drop_caches()
  raw_file = output_dir + '/' + name + '.rawdata'
  intermediate_file = output_dir + '/' + name + '.int1'
  output_file = output_dir + '/' + name + '.json'
  os.system("echo \'[\' > " + raw_file)
  print(len(binaries))
  for i in range(len(binaries) - 1):
    binary = binaries[i]
    r = 1 if (binary in bad_binaries) else rounds
    run_cmd(binary, r, graph, raw_file)

  run_cmd(binaries[len(binaries) -1], rounds, graph, raw_file)
  os.system("echo \']\' >> " + raw_file)

  os.system("cat " + raw_file  + " | grep \"^[^#;]\" > " + intermediate_file)

  os.system("""sed ':begin;N;$!bbegin;s/}\\n{/},\\n{/g' """ + intermediate_file + ' > ' + output_file)

  os.system("rm " + intermediate_file)




for graph in c_graphs:
  name = os.path.basename(graph)
  drop_caches()
  raw_file = output_dir + '/' + name + '.rawdata'
  intermediate_file = output_dir + '/' + name + '.int1'
  output_file = output_dir + '/' + name + '.json'
  os.system("echo \'[\' > " + raw_file)
  print(len(binaries))
  run_c_cmd(binaries[0], 1, graph, intermediate_file)
  for i in range(len(binaries) - 1):
    binary = binaries[i]
    r = 1 if (binary in bad_binaries) else rounds
    run_c_cmd(binary, r, graph, raw_file)

  run_c_cmd(binaries[len(binaries) -1], rounds, graph, raw_file)
  os.system("echo \']\' >> " + raw_file)

  os.system("cat " + raw_file  + " | grep \"^[^#;]\" > " + intermediate_file)

  os.system("""sed ':begin;N;$!bbegin;s/}\\n{/},\\n{/g' """ + intermediate_file + ' > ' + output_file)

  os.system("rm " + intermediate_file)
