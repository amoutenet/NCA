#pragma once

using namespace triqs::gfs;
using namespace triqs::arrays;
using mview = matrix_view<std::complex<double>>;

void solver::initialize_hybridization_functions() {

  // Hybridization functions
  Delta_gtr = triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::retime>>(gf_mesh<cartesian_product<retime,retime>>
              {{0, params.t_max, params.n_t}, {0, params.t_max, params.n_t}}, params.gf_struct);
  Delta_les = triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::retime>>(gf_mesh<cartesian_product<retime,retime>>
              {{0, params.t_max, params.n_t}, {0, params.t_max, params.n_t}}, params.gf_struct);

  // Green functions
  G_gtr = triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::retime>>(gf_mesh<cartesian_product<retime,retime>>
              {{0, params.t_max, params.n_t}, {0, params.t_max, params.n_t}}, params.gf_struct);
  G_les = triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::retime>>(gf_mesh<cartesian_product<retime,retime>>
              {{0, params.t_max, params.n_t}, {0, params.t_max, params.n_t}}, params.gf_struct);

}
