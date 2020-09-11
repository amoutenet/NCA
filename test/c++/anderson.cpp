#include <triqs/test_tools/gfs.hpp>

#include <nca/solver.hpp>
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>

using triqs::operators::c;
using triqs::operators::c_dag;
using triqs::operators::n;
using triqs::operators::many_body_operator;
using triqs::hilbert_space::gf_struct_t;

TEST(NCA, Anderson) {

  // Initialize mpi
  triqs::mpi::communicator world;

  // Parameters
  double U = 4.0;
  double beta = 2.;
  double mu = 2.0;
  double t= 1.0;
  double t_max = 2.0;
  int n_t = 10;

  // GF structure
  gf_struct_t gf_struct{{"up", {0,1}}, {"down", {0,1}}};


  // Construct CTQMC solver
  solver nca_solver({gf_struct, t_max, n_t});

  // Define complex i
  dcomplex _j {0,1};

  // Define and fill real-time Delta
  triqs::gfs::gf<refreq> r_w{{-20,20, 1000}, {1,1}};
  for (auto &w : r_w.mesh()) {
    auto omega = real(dcomplex(w));
    if (std::abs(omega) < 2) r_w[w] = 0.5 * omega - _j * std::sqrt(1 - omega * omega / 4);
    else r_w[w] = 0.5 * omega - (omega > 0 ? 1 : -1) * std::sqrt(omega * omega / 4 - 1);
  }

  triqs::gfs::gf<refreq> k_w{{-20,20, 1000}, {1,1}};
  for (auto const &w : k_w.mesh()) {
    k_w[w] = std::tanh(beta * real(dcomplex(w)) / 2.) * (r_w[w] - conj(r_w[w]));
  }

  triqs::gfs::gf<refreq> d_les_w{{-20,20, 1000}, {1,1}};
  triqs::gfs::gf<refreq> d_gtr_w{{-20,20, 1000}, {1,1}};

  for (auto const &w : r_w.mesh()) {
    auto fac = 2_j * imag(r_w[w]);
    d_les_w[w] = (k_w[w] - fac) / 2;
    d_gtr_w[w] = (k_w[w] + fac) / 2;
  }

  auto d_les_t = make_gf_from_fourier(d_les_w);
  auto d_gtr_t = make_gf_from_fourier(d_gtr_w);

  for (auto &delta: nca_solver.Delta_gtr) {
    for (auto const & t: delta.mesh()) {
	for (int a =0; a < 2; a++) {
	  delta[t](a,a) = d_gtr_t(real(dcomplex(t)))(0,0);
	}
    }
  }	

  for (auto &delta: nca_solver.Delta_les) {
    for (auto const & t: delta.mesh()) {
	for (int a =0; a < 2; a++) {
	  delta[t](a,a) = d_les_t(real(dcomplex(t)))(0,0);
	}
    }
  }	

  if (world.rank() == 0) {
    triqs::h5::file G_file("anderson.out.h5", 'w');

    h5_write(G_file, "Delta_les", nca_solver.Delta_les);
    h5_write(G_file, "Delta_gtr", nca_solver.Delta_gtr);
  }

  triqs::gfs::block_gf<triqs::gfs::retime> Delta_les, Delta_gtr;

  if (world.rank() == 0) {
    triqs::h5::file G_file("anderson.ref.h5", 'r');

    h5_read(G_file, "Delta_les", Delta_les);
    EXPECT_BLOCK_GF_NEAR(nca_solver.Delta_les, Delta_les);
    h5_read(G_file, "Delta_gtr", Delta_gtr);
    EXPECT_BLOCK_GF_NEAR(nca_solver.Delta_gtr, Delta_gtr);
  }

  // --------------------------------------------------
  
  // Define a Hamiltonian
  many_body_operator H;

  for (int i=0; i<2; i++) {
    H += U * n("up",i) * n("down",i);
    for (auto &s: {"up", "down"}) {
      H -= mu * n(s,i);
    }
  }
  for (auto &s: {"up", "down"})
    H -= t * (c_dag(s,0)*c(s,1) + c_dag(s,1)*c(s,0));

  nca_solver.initialize_atom_diag(H);

  auto R_init = nca_solver.R_gtr;
  for (int Gamma=0; Gamma<nca_solver.n_blocks; Gamma++){
      int n = nca_solver.block_sizes[Gamma];
      auto id = make_unit_matrix<std::complex<double>>(n);
      
      for (int it=n_t-1; it<2*n_t-1; it++) R_init[Gamma][it] = 1.0 * id;
  }

  nca_solver.solve({H,R_init});

  auto Z = nca_solver.get_Z();

  if (world.rank() == 0) {
    triqs::h5::file G_file("anderson.out.h5", 'a');

    h5_write(G_file, "hamilt", nca_solver.hamilt);
    h5_write(G_file, "G_les", nca_solver.G_les);
    h5_write(G_file, "G_gtr", nca_solver.G_gtr);

    h5_write(G_file, "Z", Z);
  }

  std::vector<triqs::arrays::matrix<std::complex<double>>> hamilt;
  triqs::gfs::block_gf<triqs::gfs::retime> G_gtr, G_les;
  std::complex<double> Z_check;

  if (world.rank() == 0) {
    triqs::h5::file G_file("anderson.ref.h5", 'r');
    h5_read(G_file, "hamilt", hamilt);
    h5_read(G_file, "G_les", G_les);
    h5_read(G_file, "G_gtr", G_gtr);
    h5_read(G_file, "Z", Z_check);

    for (int i=0; i<nca_solver.n_blocks; i++){
      EXPECT_ARRAY_NEAR(nca_solver.hamilt[i], hamilt[i]);
    }

    EXPECT_BLOCK_GF_NEAR(nca_solver.G_les, G_les);
    EXPECT_BLOCK_GF_NEAR(nca_solver.G_gtr, G_gtr);

    EXPECT_COMPLEX_NEAR(Z, Z_check);
  }
}

MAKE_MAIN;
