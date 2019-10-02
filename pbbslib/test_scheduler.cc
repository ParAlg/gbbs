#include "get_time.h"
#include "parse_command_line.h"
#include "utilities.h"

long fib(long i) {
  if (i <= 1) return 1;
  else if (i < 18) return fib(i-1) + fib(i-2);
  long l,r;
  par_do([&] () { l = fib(i-1);},
	 [&] () { r = fib(i-2);});
  return l + r;
}

int main (int argc, char *argv[]) {
  commandLine P(argc, argv, "[-n <size>] [-p <threads>]");
  size_t n = P.getOptionLongValue("-n", 45);
  size_t m = P.getOptionLongValue("-m", 100000000);
  size_t p = P.getOptionLongValue("-p", 0);

  auto job = [&] () {
    timer t;
    long r = fib(n);
    t.next("fib");
    cout << "result: " << r << endl;
  };

#if defined(OPENMP)
#pragma omp parallel
#pragma omp single
#endif
  parallel_run(job,p);

  auto job2 = [&] () {
    long* a = new long[m];
    auto ident = [&] (int i) {a[i] = i;};
    parallel_for(0,m,ident);
    timer t2;
    for (int i=0; i < 100; i++) {
      parallel_for(0,m,ident);
    }
    t2.next("tabulate");
  };
  parallel_run(job2,p);

  auto spin = [&] (int i) {
    for (volatile int j=0; j < 1000; j++);
  };

  auto job3 = [&] () {
    parallel_for(0,m/200,spin);
    timer t2;
    for (int i=0; i < 100; i++) {
      parallel_for(0,m/200,spin);
    }
    t2.next("map spin");
  };
  parallel_run(job3,p);
}
  
  

