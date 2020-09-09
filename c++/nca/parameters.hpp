#pragma once

#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>

// parameters of the simulation
struct constructor_p {

  // block structure of the gf
  triqs::hilbert_space::gf_struct_t gf_struct;

  // inverse temperature
  double beta;

  // number of imaginary time points
  int n_tau;

  // max real time
  double t_max;

  // numer of real time points
  int n_t;

};


struct solve_p {

  // hamiltonian
  // type: Operator
  std::function<triqs::operators::many_body_operator_generic<std::complex<double>>(double)> hamilt;

  // A boolean that tells if we have to recompute the atom diag subspaces
  bool recompute_subspaces = true;

  // convergence tolerence (not used yet)
  double tolerance = 1e-6;

  bool eq_solver = false;

};


// all parameters gathered
struct parameters: constructor_p, solve_p {

  parameters() = default;
  parameters(constructor_p const & c): constructor_p(c) {}
  parameters(constructor_p const & c, solve_p const & s): constructor_p(c), solve_p(s) {}

  void update(solve_p const & s) {
    tolerance = s.tolerance;
    hamilt = s.hamilt;
    eq_solver = s.eq_solver;
    recompute_subspaces = s.recompute_subspaces;
  }

};
