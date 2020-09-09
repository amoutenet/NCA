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
  double mu = 2.0;
  double t = 1.0;
  double beta = 2.0;
  int n_tau = 51;
  double t_max = 2.0;
  int n_t = 10;

  // GF structure
  gf_struct_t gf_struct{{"up", {0,1}}, {"down", {0,1}}};


  // Construct CTQMC solver
  solver nca_solver({gf_struct, beta, n_tau, t_max, n_t});

  // Define complex i
  dcomplex _j {0,1};

  // Define and fill imaginary-time Delta
  triqs::gfs::gf<imfreq> d_iw{{beta, Fermion, (n_tau - 1)/2},{2,2}};
  for (auto iw : d_iw.mesh()) {
    for (int a = 0; a < 2; a++) {
      auto w = imag(dcomplex(iw));
      if (w >= 0) d_iw[iw](a,a) = 0.5 * _j * w - _j * std::sqrt(1 + w*w/4);
      else d_iw[iw](a,a) = 0.5 * _j * w + _j * std::sqrt(1 + w*w/4);
    }
  }

  for (auto &delta: nca_solver.Delta_mat) {
    delta() = fourier(d_iw);
  }

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
    for (auto const & [t,tp]: delta.mesh()) {
	for (int a =0; a < 2; a++) {
	  delta[t,tp](a,a) = d_gtr_t(real(dcomplex(t)) - real(dcomplex(tp)))(0,0);
	}
    }
  }	

  for (auto &delta: nca_solver.Delta_les) {
    for (auto const & [t,tp]: delta.mesh()) {
	for (int a =0; a < 2; a++) {
	  delta[t,tp](a,a) = d_les_t(real(dcomplex(t)) - real(dcomplex(tp)))(0,0);
	}
    }
  }	

  auto r_t = make_gf_from_fourier(r_w);
  // Define and fill mixed Delta
  for (auto & delta: nca_solver.Delta_mix) {
    for (auto const & [t, tau]: delta.mesh()) {
      for (int a=0; a<2; a++) {
	delta[t, tau](a,a) = (-1) * _j * r_t(real(dcomplex(t)))(0,0) * nca_solver.Delta_mat[0](beta - real(dcomplex(tau)))(a,a);
      }
    }
  }

  // ---------------------------------------------------------------------------
  // ----------------- complete_delta() ----------------------------------------
  // ---------------------------------------------------------------------------
  
  nca_solver.complete_delta();

  if (world.rank() == 0) {
    triqs::h5::file G_file("anderson.out.h5", 'w');

    h5_write(G_file, "Delta_mat", nca_solver.Delta_mat);
    h5_write(G_file, "Delta_les", nca_solver.Delta_les);
    h5_write(G_file, "Delta_gtr", nca_solver.Delta_gtr);
    h5_write(G_file, "Delta_mix", nca_solver.Delta_mix);

    h5_write(G_file, "Delta_mix2", nca_solver.Delta_mix2);
    h5_write(G_file, "Delta_mat1", nca_solver.Delta_mat1);
    h5_write(G_file, "Delta_mat2", nca_solver.Delta_mat2);
  }

  triqs::gfs::block_gf<triqs::gfs::imtime> Delta_mat;
  triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::retime>> Delta_les, Delta_gtr;
  triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::imtime>> Delta_mix;

  triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::imtime, triqs::gfs::retime>> Delta_mix2;
  triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::imtime, triqs::gfs::imtime>> Delta_mat1, Delta_mat2;

  if (world.rank() == 0) {
    triqs::h5::file G_file("anderson.ref.h5", 'r');

    h5_read(G_file, "Delta_mat", Delta_mat);
    EXPECT_BLOCK_GF_NEAR(nca_solver.Delta_mat, Delta_mat);
    h5_read(G_file, "Delta_les", Delta_les);
    EXPECT_BLOCK_GF_NEAR(nca_solver.Delta_les, Delta_les);
    h5_read(G_file, "Delta_gtr", Delta_gtr);
    EXPECT_BLOCK_GF_NEAR(nca_solver.Delta_gtr, Delta_gtr);
    h5_read(G_file, "Delta_mix", Delta_mix);
    EXPECT_BLOCK_GF_NEAR(nca_solver.Delta_mix, Delta_mix);

    h5_read(G_file, "Delta_mix2", Delta_mix2);
    EXPECT_BLOCK_GF_NEAR(nca_solver.Delta_mix2, Delta_mix2);
    h5_read(G_file, "Delta_mat1", Delta_mat1);
    EXPECT_BLOCK_GF_NEAR(nca_solver.Delta_mat1, Delta_mat1);
    h5_read(G_file, "Delta_mat2", Delta_mat2);
    EXPECT_BLOCK_GF_NEAR(nca_solver.Delta_mat2, Delta_mat2);
  }

  // ---------------------------------------------------------------------------
  // ----------------------------- solve() -------------------------------------
  // ---------------------------------------------------------------------------

  auto function = [U, t, mu](double time) {
    // Define a Hamiltonian
    many_body_operator hamilt;
    for (int i=0; i<2; i++) {
      hamilt += U * std::cos(time) * n("up",i) * n("down",i);
      for (auto &s: {"up", "down"}) {
        hamilt -= mu * n(s,i);
      }
    }
    for (auto &s: {"up", "down"})
      hamilt -= t * (c_dag(s,0)*c(s,1) + c_dag(s,1)*c(s,0));

    return hamilt;
  };

  nca_solver.solve({function});

  auto Z = nca_solver.get_Z();
  
  if (world.rank() == 0) {
    triqs::h5::file G_file("anderson.out.h5", 'a');
    h5_write(G_file, "hamilt", nca_solver.hamilt);
    h5_write(G_file, "G_mat", nca_solver.G_mat);
    h5_write(G_file, "G_les", nca_solver.G_les);
    h5_write(G_file, "G_gtr", nca_solver.G_gtr);
    h5_write(G_file, "G_mix", nca_solver.G_mix);

    h5_write(G_file, "Z", Z);
  }

  triqs::gfs::block_gf<triqs::gfs::retime> hamilt;

  triqs::gfs::block_gf<triqs::gfs::imtime> G_mat;
  triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::retime>> G_gtr, G_les;
  triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::imtime>> G_mix;

  std::complex<double> Z_check;

  if (world.rank() == 0) {
    triqs::h5::file G_file("anderson.ref.h5", 'r');
    h5_read(G_file, "hamilt", hamilt);
    h5_read(G_file, "G_mat", G_mat);
    h5_read(G_file, "G_les", G_les);
    h5_read(G_file, "G_gtr", G_gtr);
    h5_read(G_file, "G_mix", G_mix);

    h5_read(G_file, "Z", Z_check);

    EXPECT_BLOCK_GF_NEAR(nca_solver.hamilt, hamilt);

    EXPECT_BLOCK_GF_NEAR(nca_solver.G_mat, G_mat);
    EXPECT_BLOCK_GF_NEAR(nca_solver.G_les, G_les);
    EXPECT_BLOCK_GF_NEAR(nca_solver.G_gtr, G_gtr);
    EXPECT_BLOCK_GF_NEAR(nca_solver.G_mix, G_mix);

    EXPECT_COMPLEX_NEAR(Z, Z_check, 1e-9);
  }

}

MAKE_MAIN;
