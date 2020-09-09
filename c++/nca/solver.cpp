#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/arrays/make_immutable_array.hpp>

#include "solver.hpp"
#include "volterra.hpp"
#include "nca_equation.hpp"
#include "hybridization.hpp"

#include "equilibrium.hpp"
#include "greater.hpp"
#include "mixed.hpp"
#include "lesser.hpp"

#include "green.hpp"

using namespace triqs::gfs;

solver::solver(constructor_p const & cparams) {

  // fill parameters
  params = parameters(cparams);

  // discretization times
  dt = params.t_max / (params.n_t-1);
  dtau = params.beta / (params.n_tau-1);

  // determine basis of operators to use from gf_struct
  for (auto const &bl : params.gf_struct) {
    for (auto const &a : bl.second) { fops.insert(bl.first, a); }
  }

  // initialize hybridization functions
  initialize_hybridization_functions();

}


void solver::solve(solve_p const & sparams) {

  // update parameters
  params.update(sparams);

  if (params.recompute_subspaces) initialize_atom_diag(params.hamilt);
  update_hamiltonian(params.hamilt);

  // print some info
  if (world.rank() == 0) {

    std::cout << "Welcome to the NCA solver" << std::endl << std::endl;
    if (params.eq_solver) std::cout << "You are only solving the equilibrium part of the contour" << std::endl;
    std::cout << "Number of subspaces: " << atom.n_subspaces() << std::endl;
    std::cout << "Subspace sizes: ";
    for (int i=0; i<n_blocks; i++) std::cout << block_sizes[i] << " ";
    std::cout << std::endl;
    std::cout << "Number of particles: ";
    for (int i=0; i<n_blocks; i++) std::cout << n_particles[i] << " ";
    std::cout << std::endl;
    std::cout << "dt = " << dt << std::endl;
    std::cout << "dtau = " << dtau << std::endl << std::endl;

  }

  // complete the hybridization functions
  complete_delta();

  // solve all components
  solve_equilibrium();
  if (!params.eq_solver) {
    solve_greater();
    solve_mixed();
    solve_lesser();
  }

  // compute Greens' functions
  compute_G_mat();
  if (!params.eq_solver) {
    compute_G_mix();
    compute_G_les();
    compute_G_gtr();
  }

}


void solver::initialize_atom_diag(std::function<triqs::operators::many_body_operator_generic<std::complex<double>>(double)> function) {
  // initialize atomic diagonalization
  atom = {function(0), fops};

  n_blocks = atom.n_subspaces();
  block_sizes.resize(n_blocks);
  n_particles.resize(n_blocks);

  // get info from diagonalization
  for (int Gamma=0; Gamma<n_blocks; Gamma++) block_sizes[Gamma] = atom.get_subspace_dim(Gamma);

  // count number of particles in every block
  for (int Gamma=0; Gamma<n_blocks; Gamma++) {
    int n = atom.get_fock_states()[Gamma][0];
    int count = 0;
    while (n) {
      count += n & 1;
      n >>= 1;
    }
    n_particles[Gamma] = count;
  }

  // initialize R and S for all components
  initialize_R_and_S();

  // initialize hamiltonian
  initialize_hamiltonian();
}


void solver::initialize_hamiltonian() {

  // Construct the vector of gf
  std::vector<gf<retime>> gf_vec;
  
  // Initialize Hamiltonian
  for (int Gamma=0; Gamma<n_blocks; Gamma++) {
    int n = block_sizes[Gamma];
    gf_vec.emplace_back(gf_mesh<retime>{0, params.t_max, params.n_t}, make_shape(n, n));
  }
  
  hamilt = block_gf<retime>(std::move(gf_vec));
  
}
  

void solver::update_hamiltonian(std::function<triqs::operators::many_body_operator_generic<std::complex<double>>(double)> function) {

  for (int i=0; i<params.n_t; i++) {
    auto op_mat = atom.get_op_mat(function(i*dt));

    for (int Gamma=0; Gamma<n_blocks; Gamma++) {
      if (op_mat.connection(Gamma) == -1) hamilt[Gamma][i] = 0;
      else hamilt[Gamma][i] = op_mat.block_mat[Gamma];
    }
  }
}


void solver::initialize_R_and_S() {

  // Construct the vector of gf for equilibrium
  std::vector<gf<imtime>> gf_vec_tau;

  for (int Gamma=0; Gamma<n_blocks; Gamma++) {
    int n = block_sizes[Gamma];
    gf_vec_tau.emplace_back(gf_mesh<imtime>{params.beta, Fermion, params.n_tau}, make_shape(n, n));
  }

  R_eq = block_gf<imtime>(gf_vec_tau);
  Rdot_eq = block_gf<imtime>(gf_vec_tau);
  S_eq = block_gf<imtime>(gf_vec_tau);


  // Construct the vector of gf for lesser and greater
  std::vector<gf<cartesian_product<retime, retime>>> gf_vec_2t;

  for (int Gamma=0; Gamma<n_blocks; Gamma++) {
    int n = block_sizes[Gamma];
    gf_vec_2t.emplace_back(gf_mesh<cartesian_product<retime, retime>>{{0, params.t_max, params.n_t},{0, params.t_max, params.n_t}}, make_shape(n, n));
  }

  R_gtr = block_gf<cartesian_product<retime, retime>>(gf_vec_2t);
  Rdot_gtr = block_gf<cartesian_product<retime, retime>>(gf_vec_2t);
  S_gtr = block_gf<cartesian_product<retime, retime>>(gf_vec_2t);

  R_les = block_gf<cartesian_product<retime, retime>>(gf_vec_2t);
  Rdot_les = block_gf<cartesian_product<retime, retime>>(gf_vec_2t);
  S_les = block_gf<cartesian_product<retime, retime>>(gf_vec_2t);
  Q_les = block_gf<cartesian_product<retime, retime>>(gf_vec_2t);


  // Construct the vector of gf for mixed
  std::vector<gf<cartesian_product<retime, imtime>>> gf_vec_ttau;

  for (int Gamma=0; Gamma<n_blocks; Gamma++) {
    int n = block_sizes[Gamma];
    gf_vec_ttau.emplace_back(gf_mesh<cartesian_product<retime, imtime>>{{0, params.t_max, params.n_t},{params.beta, Fermion, params.n_tau}}, make_shape(n, n));
  }

  R_mix = block_gf<cartesian_product<retime, imtime>>(gf_vec_ttau);
  Rdot_mix = block_gf<cartesian_product<retime, imtime>>(gf_vec_ttau);
  S_mix = block_gf<cartesian_product<retime, imtime>>(gf_vec_ttau);
  Q_mix = block_gf<cartesian_product<retime, imtime>>(gf_vec_ttau);

}

