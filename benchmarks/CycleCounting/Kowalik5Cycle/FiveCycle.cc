#include "FiveCycle.h"
#include "FiveCycle-escape.h"
#include "FiveCycle-match.h"

long strToAlg(std::string alg_str){
  if (alg_str == "KOWALIK") return 0;
  else if (alg_str == "ESCAPE") return 1;
  else if (alg_str == "EXP") return 2;
  else if (alg_str == "MATCH") return 3;
  std::cout << "Using default algorithm (Kowalik)" << std::endl;
  return 0;
}


template <class Graph>
double Count5Cycle_runner(Graph& G, commandLine P) {
  std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
  std::cout << "### Application: Cycle" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### ------------------------------------" << endl;
  assert(P.getOption("-s")); // make sure input graph is symmetric
  double epsilon = P.getOptionDoubleValue("-e", 0.1); // epsilon, for Goodrich-Pszona or Barenboim-Elkin directing
  bool serial = P.getOptionValue("--serial"); 
  auto alg_str = P.getOptionValue("--alg", "");
  long alg = strToAlg(alg_str);
  if (alg == 0) alg_str = "KOWALIK"; 

  bool no_schedule = P.getOptionValue("--no-schedule");
  // bool experiment = P.getOptionValue("--exp");
  // bool escape = P.getOptionValue("--escape");
  long order_type = P.getOptionLongValue("-o", 0); 
  long block_size = P.getOptionLongValue("-b", 500000);
  std::string dirtype_str;
  if (order_type == 0){
    dirtype_str = "Goodrich-Pszona";
  } else if (order_type == 1) {
    dirtype_str = "Barenboim-Elkin";
  } else if (order_type == 2) {
    dirtype_str = "Degree";
  } else if (order_type == 3) {
    dirtype_str = "KCore";
  }
  std::cout << "### Direct Type: " << dirtype_str << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  ulong numCycles;
  long block_size_display = (alg == 2) || serial || no_schedule? -1 : block_size; 
  std::cout << "##### Experiment: alg:"<< alg_str << ",par:" << !serial << ",threads:" << num_workers() << ",sched:" << (!no_schedule & !serial) << ",bs:" << block_size_display << "," << dirtype_str << "," << P.getArgument(0) << std::endl;
  timer t; t.start();
  if (alg == 3) {
    numCycles = Count5Cycle_match(G, order_type, block_size, epsilon);
  } else if (alg == 2) { // experimental
    numCycles = Count5Cycle_experiment(G, order_type, block_size, epsilon);
  } else if (alg == 1) { // ESCAPE
    if (serial){
      numCycles = Count5Cycle_ESCAPE(G, order_type, block_size, epsilon);
    } else if (no_schedule) {
      // block_size_display = block_size;
      numCycles = Count5Cycle_ESCAPE_no_scheduling(G, order_type, block_size, epsilon);
    } else {
      block_size_display = block_size;
      numCycles = Count5Cycle_ESCAPE_par(G, order_type, block_size, epsilon);
    }
  } else { // Kowalik
    if (serial){
      numCycles = Count5Cycle_serial(G, order_type, block_size, epsilon);
    } else if (no_schedule) { 
      numCycles = Count5Cycle_no_scheduling(G, order_type, block_size, epsilon);
    }else {
      block_size_display = block_size;
      numCycles = Count5Cycle(G, order_type, block_size, epsilon);
    }
  }
  double tt = t.stop();
  
  std::cout << "### Number of Cycles: " << numCycles << std::endl;
  std::cout << "### Running Time: " << tt << std::endl;
  return tt;
}

generate_symmetric_main(Count5Cycle_runner, false);