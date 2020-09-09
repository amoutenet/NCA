#pragma once

using namespace triqs::arrays;

void solver::solve_mixed() {

  auto R_mix_old = R_mix;

  // first column t = 0 for R_mix
  for (int Gamma=0; Gamma<n_blocks; Gamma++) {
    for (int tau=0; tau<params.n_tau; tau++)
      R_mix[Gamma][{0,tau}] = 1_j * R_eq[Gamma][{params.n_tau-1-tau}];
  }

  // first column t = 0 for Rdot_mix
  for (int Gamma=0; Gamma<n_blocks; Gamma++) {

    for (int tau=0; tau<params.n_tau; tau++)
      get_S_from_NCA(Gamma, 0, tau, S_mix, R_mix, Delta_mix, Delta_mix2, -1);

    for (int tau=0; tau<params.n_tau; tau++) {
      get_mixed_Q(Gamma, 0, tau);
      Rdot_mix[Gamma][{0,tau}] = -1_j * hamilt[Gamma][0] * R_mix[Gamma][{0,tau}] - 1_j * Q_mix[Gamma][{0,tau}];
    }

  }

  // move in time
  int nca_loops = 0;
  for (int it=1; it<params.n_t; it++) {

    // initial guess for R
    for (int Gamma=0; Gamma<n_blocks; Gamma++){
        for (int tau=0; tau<params.n_tau; tau++) {
          R_mix[Gamma][{it,tau}] = 1.0;
	}
    }

    R_mix_old = R_mix;

    // NCA loop
    std::vector<double> conv(n_blocks);
    conv[0] = 1;

    nca_loops = 0;
    while (*std::max_element(conv.begin(), conv.end()) > params.tolerance) {
      nca_loops += 1;

      // get R by solving Volterra equation for all blocks
      for (int Gamma=0; Gamma<n_blocks; Gamma++) {

        // get S from NCA equations and then Q
        for (int tau=0; tau<params.n_tau; tau++)
          get_S_from_NCA(Gamma, it, tau, S_mix, R_mix, Delta_mix, Delta_mix2, -1);

        for (int tau=0; tau<params.n_tau; tau++) {
          get_mixed_Q(Gamma, it, tau);
          volterra_step(R_mix[Gamma], Rdot_mix[Gamma], S_gtr[Gamma], hamilt[Gamma], Q_mix[Gamma], dt, it, tau, 0);
        }

      }

      for (int Gamma=0; Gamma<n_blocks; Gamma++) {
        conv[Gamma] = max_element(abs(R_mix[Gamma].data() - R_mix_old[Gamma].data())) / max_element(abs(R_mix[Gamma].data() + R_mix_old[Gamma].data()));
      }

      R_mix_old = R_mix;
    } // nca loop

  } // time loop

  std::cout << "Mixed part finished" << std::endl;
  std::cout << "Numver of NCA loops for last time point: " << nca_loops << std::endl << std::endl;

}


// helper function
void solver::get_mixed_Q(int Gamma, int t, int tau) {

  Q_mix[Gamma][{t,tau}] = 0.0;
  for (int k=tau; k<params.n_tau-1; k++) {
    Q_mix[Gamma][{t,tau}] += 0.5 * dtau * (
                 S_mix[Gamma][{t,k}] * R_eq[Gamma][k-tau] +
                 S_mix[Gamma][{t,k+1}] * R_eq[Gamma][k+1-tau] );
  }

}
