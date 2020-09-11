#pragma once

#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>

// parameters of the simulation
struct constructor_p {

  // block structure of the gf
  triqs::hilbert_space::gf_struct_t gf_struct;

  // max real time
  double t_max;

  // numer of real time points
  int n_t;

};


struct solve_p {

  std::vector<triqs::arrays::array<std::complex<double>,2>> R_init;

  // hamiltonian
  // type: Operator
  std::function<triqs::operators::many_body_operator_generic<std::complex<double>>(double)> hamilt;

  bool recompute_subspaces = true;

  // convergence tolerence (not used yet)
  double tolerance = 1e-6;

};


// all parameters gathered
struct parameters: constructor_p, solve_p {

  parameters() = default;
  parameters(constructor_p const & c): constructor_p(c) {}
  parameters(constructor_p const & c, solve_p const & s): constructor_p(c), solve_p(s) {}

  void update(solve_p const & s) {
    tolerance = s.tolerance;
    R_init = s.R_init;
    hamilt = s.hamilt;
    recompute_subspaces = s.recompute_subspaces;
  }

};
