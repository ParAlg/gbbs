#pragma once
#include "basic_node.h"

// *******************************************
//   AUGMENTED NODE
// *******************************************

// Creates an augmented node from a basic node.
// The new augmented entry is a pair of the original entry and the aumented
// value.   The Entry struct must have the static inteface:
//   entry_t;
//   aug_t;
//   get_empty() -> aug_t;
//   from_entry(entry_t) -> aug_t;
//   combine(aut_t, aug_t) -> aug_t;
template<class balance, class Entry>
struct aug_node : basic_node<balance, std::pair<typename Entry::entry_t,
						typename Entry::aug_t>> {
  using AT = typename Entry::aug_t;
  using ET = typename Entry::entry_t;
  using basic = basic_node<balance, std::pair<ET,AT>>;
  using node = typename basic::node;

  static ET& get_entry(node *a) {return a->entry.first;}
  static ET* get_entry_p(node *a) {return &a->entry.first;}
  static void set_entry(node *a, ET e) {a->entry.first = e;}

  static AT aug_val(node* a) {
    if (a == NULL) return Entry::get_empty();
    else return (a->entry).second;}

  static void update(node* a) {
    basic::update(a);
    AT av = Entry::from_entry(get_entry(a));
    if (a->lc) av = Entry::combine(((a->lc)->entry).second, av);
    if (a->rc) av = Entry::combine(av, ((a->rc)->entry).second);
    (a->entry).second = av;
  }

  // updates augmented value using f, instead of recomputing
  template<class F>
  static void lazy_update(node* a, F f) {
    basic::update(a);
    (a->entry).second = f((a->entry).second);
  }

  static node* make_node(ET e) {
    std::pair<ET,AT> ea;
    ea.first = e;
    return basic::make_node(ea);
  }

  static node* single(ET e) {
    AT av = Entry::from_entry(e);
    return basic::single(std::make_pair(e,av));
  }
};
