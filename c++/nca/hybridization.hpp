#pragma once

using namespace triqs::gfs;

void solver::initialize_hybridization_functions() {

  // Hybridization functions
  Delta_gtr = block_gf<cartesian_product<retime,retime>>(gf_mesh<cartesian_product<retime,retime>>
              {{0, params.t_max, params.n_t}, {0, params.t_max, params.n_t}}, params.gf_struct);
  Delta_les = block_gf<cartesian_product<retime,retime>>(gf_mesh<cartesian_product<retime,retime>>
              {{0, params.t_max, params.n_t}, {0, params.t_max, params.n_t}}, params.gf_struct);

  Delta_mix = block_gf<cartesian_product<retime,imtime>>(gf_mesh<cartesian_product<retime,imtime>>
              {{0, params.t_max, params.n_t}, {params.beta, Fermion, params.n_tau}}, params.gf_struct);
  Delta_mix2 = block_gf<cartesian_product<imtime,retime>>(gf_mesh<cartesian_product<imtime,retime>>
              {{params.beta, Fermion, params.n_tau}, {0, params.t_max, params.n_t}}, params.gf_struct);

  Delta_mat = block_gf<imtime>(gf_mesh<imtime>{params.beta, Fermion, params.n_tau}, params.gf_struct);

  Delta_mat1 = block_gf<cartesian_product<imtime,imtime>>(gf_mesh<cartesian_product<imtime,imtime>>
              {{params.beta, Fermion, params.n_tau}, {params.beta, Fermion, 1}}, params.gf_struct);
  Delta_mat2 = block_gf<cartesian_product<imtime,imtime>>(gf_mesh<cartesian_product<imtime,imtime>>
              {{params.beta, Fermion, 1}, {params.beta, Fermion, params.n_tau}}, params.gf_struct);

  // Green functions
  G_mat = block_gf(gf_mesh<imtime>{params.beta, Fermion, params.n_tau}, params.gf_struct);

  G_gtr = block_gf<cartesian_product<retime,retime>>(gf_mesh<cartesian_product<retime,retime>>
          {{0, params.t_max, params.n_t}, {0, params.t_max, params.n_t}}, params.gf_struct);
  G_les = block_gf<cartesian_product<retime,retime>>(gf_mesh<cartesian_product<retime,retime>>
          {{0, params.t_max, params.n_t}, {0, params.t_max, params.n_t}}, params.gf_struct);

  G_mix = block_gf<cartesian_product<retime,imtime>>(gf_mesh<cartesian_product<retime,imtime>>
          {{0, params.t_max, params.n_t}, {params.beta, Fermion, params.n_tau}}, params.gf_struct);

}


// The functions below are necessary to use the generic get_S_from_NCA
void solver::complete_delta() {

  int i = 0;
  for (auto const &bl : params.gf_struct) {

    for (int tau=0; tau<params.n_tau; tau++)
      for (int t=0; t<params.n_t; t++)
        Delta_mix2[i][{params.n_tau-1-tau,t}] = dagger(triqs::make_view(Delta_mix[i][{t,tau}]));

    for (int tau=0; tau<params.n_tau; tau++) {
      Delta_mat1[i][{tau,0}] = Delta_mat[i][tau];
      Delta_mat2[i][{0,tau}] = -Delta_mat[i][params.n_tau-1-tau];
    }

    i++;

  }

}

