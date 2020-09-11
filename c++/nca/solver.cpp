#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/arrays/make_immutable_array.hpp>

#include "solver.hpp"
#include "nca_equation.hpp"
#include "hybridization.hpp"

#include "greater.hpp"
#include "lesser.hpp"

#include "green.hpp"

solver::solver(constructor_p const & cparams) {

  // fill parameters
  params = parameters(cparams);

  // discretization times
  dt = params.t_max / (params.n_t-1);

  // determine basis of operators to use from gf_struct
  for (auto const & bl : params.gf_struct) {
    for (auto const & a : bl.second) {fops.insert(bl.first, a);}
  }
  // Why is the command below not compiling??
  //for (auto const & [bname, idx_lst] : params.gf_struct) fops.insert(bname, idx_lst);

  // initialize hybridization functions
  initialize_hybridization_functions();

  // initialize vector needed to treat singularities
  rho = std::vector<double>(params.n_t);

}


void solver::solve(solve_p const & sparams) {
  // update parameters
  params.update(sparams);

  if (params.recompute_subspaces) initialize_atom_diag(params.H);
  update_hamiltonian(params.H);

  // print some info
  if (world.rank() == 0) {

    std::cout << "Welcome to the NESS NCA solver" << std::endl << std::endl;
    std::cout << "Number of subspaces: " << atom.n_subspaces() << std::endl;
    std::cout << "Subspace sizes: ";
    for (int i=0; i<n_blocks; i++) std::cout << block_sizes[i] << " ";
    std::cout << std::endl;
    std::cout << "Number of particles: ";
    for (int i=0; i<n_blocks; i++) std::cout << n_particles[i] << " ";
    std::cout << std::endl;
    std::cout << "dt = " << dt << std::endl;

  }


  // solve all components
  solve_greater();
  solve_lesser();

  compute_G_gtr();
  compute_G_les();
}


void solver::initialize_atom_diag(triqs::operators::many_body_operator_generic<std::complex<double>> H) {
  // initialize atomic diagonalization
  atom = {H, fops};

  // get info from diagonalization
  n_blocks = atom.n_subspaces();
  block_sizes.resize(n_blocks);
  n_particles.resize(n_blocks);

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


void solver::update_hamiltonian(triqs::operators::many_body_operator_generic<std::complex<double>> H) {

  auto op_mat = atom.get_op_mat(H);

  for (int Gamma=0; Gamma<n_blocks; Gamma++) {
      if (op_mat.connection(Gamma) == -1) hamilt[Gamma]() = 0;
      else hamilt[Gamma] = op_mat.block_mat[Gamma];
  }
}

void solver::initialize_hamiltonian() {

  hamilt.resize(0);

  // Initialize Hamiltonian
  for (int Gamma=0; Gamma<n_blocks; Gamma++) {
    int n = block_sizes[Gamma];
    hamilt.emplace_back(n,n);
  }

}


void solver::initialize_R_and_S() {

  // Construct the vector of gf
  std::vector<triqs::gfs::gf<triqs::gfs::retime>> gf_vec;

  for (int Gamma=0; Gamma<n_blocks; Gamma++) {
    int n = block_sizes[Gamma];
    gf_vec.emplace_back(triqs::gfs::gf_mesh<triqs::gfs::retime>{-params.t_max, params.t_max, 2*params.n_t - 1}, triqs::gfs::make_shape(n, n));
  }

  // Initialize
  R_gtr = triqs::gfs::block_gf<triqs::gfs::retime>(gf_vec);
  R_tilde_gtr = triqs::gfs::block_gf<triqs::gfs::retime>(gf_vec);
  S_gtr = triqs::gfs::block_gf<triqs::gfs::retime>(gf_vec);
  S_tilde_gtr = triqs::gfs::block_gf<triqs::gfs::retime>(gf_vec);

  R_les = triqs::gfs::block_gf<triqs::gfs::retime>(gf_vec);
  S_les = triqs::gfs::block_gf<triqs::gfs::retime>(gf_vec);

  fit_R = triqs::gfs::block_gf<triqs::gfs::retime>(gf_vec);
  temp_R = triqs::gfs::block_gf<triqs::gfs::retime>(gf_vec);
  fit_S = triqs::gfs::block_gf<triqs::gfs::retime>(gf_vec);
  temp_S = triqs::gfs::block_gf<triqs::gfs::retime>(gf_vec);

  R_tilde_gtr_w = make_gf_from_fourier(R_tilde_gtr);
  S_tilde_gtr_w = make_gf_from_fourier(S_tilde_gtr);
  R_gtr_w = make_gf_from_fourier(R_gtr);
  S_gtr_w = make_gf_from_fourier(S_gtr);
  R_les_w = make_gf_from_fourier(R_les);
  S_les_w = make_gf_from_fourier(S_les);

  fit_R_w = make_gf_from_fourier(fit_R);
  temp_R_w = make_gf_from_fourier(temp_R);
  fit_S_w = make_gf_from_fourier(fit_S);
  temp_S_w = make_gf_from_fourier(temp_S);

}

