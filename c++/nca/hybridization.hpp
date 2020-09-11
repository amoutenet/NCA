#pragma once

using namespace triqs::gfs;

void solver::initialize_hybridization_functions() {

  // Hybridization functions
  Delta_gtr = block_gf<retime>(gf_mesh<retime>{-params.t_max, params.t_max, 2*params.n_t - 1}, params.gf_struct);
  Delta_les = block_gf<retime>(gf_mesh<retime>{-params.t_max, params.t_max, 2*params.n_t - 1}, params.gf_struct);

  // Green functions
  G_les = block_gf<retime>(gf_mesh<retime>{-params.t_max, params.t_max, 2*params.n_t - 1}, params.gf_struct);
  G_gtr = block_gf<retime>(gf_mesh<retime>{-params.t_max, params.t_max, 2*params.n_t - 1}, params.gf_struct);

}
