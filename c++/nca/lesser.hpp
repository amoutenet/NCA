#pragma once

using namespace triqs::arrays;

void solver::solve_lesser() {

  // Initial guess for R
  for (int Gamma=0; Gamma<n_blocks; Gamma++){
    int n = block_sizes[Gamma];
    auto id = make_unit_matrix<std::complex<double>>(n);

    for (int it=params.n_t-1; it<2*params.n_t-1; it++) R_les[Gamma][it] = 1.0 * id;
  }

  auto R_les_old = R_les;

  // NCA loop
  std::vector<double> conv_abs(n_blocks);
  std::vector<double> conv_rel(n_blocks);
  double err_rel = 1.;

  int nca_loops = 0;
  while ((err_rel > params.tolerance) && (nca_loops < 200)) {
    nca_loops += 1;

    /*
    std::cout << "NCA loop number " << nca_loops << std::endl;
    for (int i=0; i<conv_abs.size(); i++) {
      std::cout << conv_abs[i] << " - ";
    }
    std::cout << std::endl;
    for (int i=0; i<conv_rel.size(); i++) {
      std::cout << conv_rel[i] << " - ";
    }
    std::cout << std::endl;
    */

    std::complex<double> Z = 0;
    for (int Gamma=0; Gamma<n_blocks; Gamma++) {
      int n = block_sizes[Gamma];
      auto id = make_unit_matrix<std::complex<double>>(n);

      auto km = make_zero_tail(R_les_w[Gamma], 4);

      // get S_les[Gamma] from R_les
      get_S_from_NCA(Gamma, S_les, R_les, Delta_les, Delta_gtr, -1);

      // fill negative time points
      for (int it = 0; it<params.n_t - 1; it++){
	S_les[Gamma][it] = -dagger(triqs::make_view(S_les[Gamma][2*params.n_t-2-it]));
      }

      S_les_w[Gamma] = make_gf_from_fourier(S_les[Gamma], S_les_w[Gamma].mesh(), km);

      // get R_les[Gamma] from S_les
      for (auto const & w : R_les_w[Gamma].mesh()) {
        R_les_w[Gamma][w] = (-1_j) * imag(inverse(real(dcomplex(w))* id - hamilt[Gamma] - S_tilde_gtr_w[Gamma][w]) * S_les_w[Gamma][w] * R_tilde_gtr_w[Gamma][w]);
      }

      R_les[Gamma] = make_gf_from_fourier(R_les_w[Gamma], R_les[Gamma].mesh(), km);
      Z += 1_j * trace(R_les[Gamma][params.n_t-1]);
    }

    for (int Gamma = 0; Gamma < n_blocks; Gamma++){
      R_les[Gamma] = R_les[Gamma] / Z;
      conv_abs[Gamma] = max_element(abs(R_les[Gamma].data() - R_les_old[Gamma].data()));

      if (conv_abs[Gamma] < 1e-9) conv_rel[Gamma] = 0.;
      else conv_rel[Gamma] = max_element(abs(R_les[Gamma].data() - R_les_old[Gamma].data())) / max_element(abs(R_les[Gamma].data() + R_les_old[Gamma].data()));
    }
    
    err_rel = *std::max_element(conv_rel.begin(), conv_rel.end());

    R_les_old = R_les;

  } // nca loop

  std::cout << "Lesser part finished" << std::endl;
  std::cout << "Number of NCA loops: " << nca_loops << std::endl << std::endl;

}
