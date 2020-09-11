#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/arrays/make_immutable_array.hpp>

#include "solver.hpp"
#include "volterra.hpp"
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

}


void solver::solve(solve_p const & sparams) {
  // update parameters
  params.update(sparams);

  if (params.recompute_subspaces) initialize_atom_diag(params.hamilt);
  update_hamiltonian(params.hamilt);

  // print some info
  if (world.rank() == 0) {

    std::cout << "Welcome to the NCA solver computing the transient on the two branch contour" << std::endl << std::endl;
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
  solve_lesser(params.R_init);

  compute_G_gtr();
  compute_G_les();
}


void solver::initialize_atom_diag(std::function<triqs::operators::many_body_operator_generic<std::complex<double>>(double)> function) {
  // initialize atomic diagonalization
  atom = {function(0), fops};

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


void solver::update_hamiltonian(std::function<triqs::operators::many_body_operator_generic<std::complex<double>>(double)> function) {

  for (int it=0; it < params.n_t; it++) {
    auto op_mat = atom.get_op_mat(function(it * dt));

    for (int Gamma=0; Gamma<n_blocks; Gamma++) {
      if (op_mat.connection(Gamma) == -1) hamilt[Gamma][it] = 0;
      else hamilt[Gamma][it] = op_mat.block_mat[Gamma];
    }
  }
}

void solver::initialize_hamiltonian() {

  // Construct the vector of gf
  std::vector<triqs::gfs::gf<triqs::gfs::retime>> gf_vec;

  // Initialize Hamiltonian
  for (int Gamma=0; Gamma<n_blocks; Gamma++) {
    int n = block_sizes[Gamma];
    gf_vec.emplace_back(triqs::gfs::gf_mesh<triqs::gfs::retime>{0, params.t_max, params.n_t}, triqs::gfs::make_shape(n, n));
  }

  hamilt = triqs::gfs::block_gf<triqs::gfs::retime>(std::move(gf_vec));

}


void solver::initialize_R_and_S() {

  // Construct the vector of gf
  std::vector<triqs::gfs::gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::retime>>> gf_vec;

  for (int Gamma=0; Gamma<n_blocks; Gamma++) {
    int n = block_sizes[Gamma];
    gf_vec.emplace_back(triqs::gfs::gf_mesh<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::retime>>{{0, params.t_max, params.n_t},{0, params.t_max, params.n_t}}, triqs::gfs::make_shape(n, n));
  }

  // Initialize
  R_gtr = triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::retime>>(gf_vec);
  Rdot_gtr = triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::retime>>(gf_vec);
  S_gtr = triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::retime>>(gf_vec);

  R_les = triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::retime>>(gf_vec);
  Rdot_les = triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::retime>>(gf_vec);
  S_les = triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::retime>>(gf_vec);
  Q_les = triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::retime>>(gf_vec);

}

