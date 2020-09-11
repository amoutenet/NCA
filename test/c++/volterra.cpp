#include <triqs/test_tools/gfs.hpp>

#include <nca/solver.hpp>
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>

using triqs::operators::c;
using triqs::operators::c_dag;
using triqs::operators::n;
using triqs::hilbert_space::gf_struct_t;
using triqs::arrays::matrix;
using triqs::arrays::range;

using namespace triqs::gfs;

// A non-trivial matrix case
TEST(NCA, VolterraMatrix) {

  // Initialize mpi
  triqs::mpi::communicator world;

  // Parameters
  double beta = 1.0;
  double t_max = 10.0;
  int n_t = 100;
  double U = 10;

  // GF structure
  gf_struct_t gf_struct{{"0", {0,1}}};
  auto hamilt = c_dag("0",1) * c("0",0) + c_dag("0",0) * c("0",1) + U * (n("0",0) + n("0",1));
  auto function = [hamilt](double time) {return hamilt;};

  int n = 2;

  auto H = gf<retime>{gf_mesh<retime>{0, t_max, n_t}, make_shape(n, n)};
  auto Q = gf<cartesian_product<retime, retime>>{gf_mesh<cartesian_product<retime, retime>>{{0, t_max, n_t},{0, t_max, n_t}}, make_shape(n, n)};

  // Construct CTQMC solver
  solver ns({gf_struct, t_max, n_t});

  ns.initialize_atom_diag(function);
     
  // initial conditions
  for (int k=0; k<n; k++) {
    for (int l=0; l<n; l++) {

        // Set non-trivial H, S and Q
        for (int i=0; i<n_t; i++) {

          double taur = i*ns.dt;
          H[i](k,l) = (0.5 * 1_j + 0.1) * cos(taur) + k * 0.2 - l * 0.4 * 1_j;
          Q[{i,0}](k,l) = (-0.04 * 1_j + 0.1) * cos(4*taur) + (0.3*k - 0.1*l* 1_j) * sin(2*taur);

          for (int j=0; j<i+1; j++) {
            double taud = (i-j)*ns.dt;
            ns.S_gtr[1][{i,j}](k,l) = -0.5 * 1_j * (exp(-taud) + cos(taud))
                                  + 0.2 * (exp(-2*taud) + sin(4*taud));
          }

        }

        // Initial conditions
        ns.R_gtr[1][{0,0}](k,l) = -1 + 0.5 * 1_j + k * 0.2 + l * 0.4 * 1_j;

    }
  }

  ns.Rdot_gtr[1][{0,0}] = -1_j * H[0] * ns.R_gtr[1][{0,0}] -1_j * Q[{0,0}];

  // solve Volterra equation
  for (int t=1; t<n_t; t++)
    ns.volterra_step(ns.R_gtr[1], ns.Rdot_gtr[1], ns.S_gtr[1], H, Q, ns.dt, t, 0, 0);


  // write results to file
  triqs::h5::file f("volterra.out.h5", 'a');
  h5_write(f, "R_2", ns.R_gtr[1]);
  h5_write(f, "Rdot_2", ns.Rdot_gtr[1]);

  // read reference
  gf<cartesian_product<retime, retime>> R_check, Rdot_check;
  triqs::h5::file fc("volterra.ref.h5", 'r');
  h5_read(fc, "R_2", R_check);
  h5_read(fc, "Rdot_2", Rdot_check);

  // compare
  EXPECT_GF_NEAR(ns.R_gtr[1], R_check);
  EXPECT_GF_NEAR(ns.Rdot_gtr[1], Rdot_check);

}

MAKE_MAIN;
