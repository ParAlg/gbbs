#pragma once

/* ************************* Benchmark Utils *************************** */
namespace gbbs {
static pbbs::timer bt;
using uchar = unsigned char;

#define time(_var,_body)    \
  bt.start();               \
  _body;		    \
  double _var = bt.stop();

namespace benchmark {
  template<typename F>
  double reduce(std::vector<double> V, F f) {
    double x = V[0];
    for (size_t i=1; i < V.size(); i++) x = f(x,V[i]);
    return x;
  }

  template <class T>
  double median(std::vector<T> V) {
    std::sort(V.begin(),V.end());
    if (V.size()%2 == 1)
      return V[V.size()/2];
    else
      return (V[V.size()/2] + V[V.size()/2 - 1])/2.0;
  }

  double sumf(double a, double b) {return a+ b;};
  double minf(double a, double b) {return (a < b) ? a : b;};
  double maxf(double a, double b) {return (a > b) ? a : b;};

  /* returns a vector of running times for test(), run rounds many times. */
  template<typename F>
  std::vector<double> repeat(size_t rounds, F& test) {
    std::vector<double> R;
    for (size_t i=0; i < rounds; i++) {
      auto t = test();
      std::cout << "### t = " << t << std::endl;
      R.push_back(t);
    }
    return R;
  }

  /* returns a triple of the (minimum, maximum, medium) running times for test */
  template<typename F>
  std::tuple<double, double, double> run_multiple(size_t rounds, F& test) {
    std::vector<double> t = repeat(rounds, test);

    double mint = reduce(t, minf);
    double maxt = reduce(t, maxf);
    double medt = median(t);
    return std::make_tuple(mint, maxt, medt);
  }
}
}  // namespace gbbs
