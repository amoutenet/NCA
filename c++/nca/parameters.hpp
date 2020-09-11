#pragma once

#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>

// parameters of the simulation
struct constructor_p {

  // block structure of the gf
  triqs::hilbert_space::gf_struct_t gf_struct;

  // max real time
  double t_max;

  // numer of positive real time points (total number of pts: 2*n_t)
  int n_t;

};


struct solve_p {

  // hamiltonian
  // type: Operator
  triqs::operators::many_body_operator_generic<std::complex<double>> H;

  // First guess for R^>
  triqs::gfs::block_gf<triqs::gfs::retime> R_init;

  // do we want to recompute the subspace from atom_diag class?
  bool recompute_subspaces = true;

  // convergence tolerence
  double tolerance = 1e-6;

  // threshold to determine the support of S
  double S_threshold = 1e-9;
};


// all parameters gathered
struct parameters: constructor_p, solve_p {

  parameters() = default;
  parameters(constructor_p const & c): constructor_p(c) {}
  parameters(constructor_p const & c, solve_p const & s): constructor_p(c), solve_p(s) {}

  void update(solve_p const & s) {
    H = s.H;
    R_init = s.R_init;
    recompute_subspaces = s.recompute_subspaces;
    tolerance = s.tolerance;
    S_threshold = s.S_threshold;

  }

};
