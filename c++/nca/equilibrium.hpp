#pragma once

using namespace triqs::arrays;

void solver::solve_equilibrium() {

  // initial condition and guess
  auto R_eq_old = R_eq;

  for (int Gamma=0; Gamma<n_blocks; Gamma++) {

    int n = block_sizes[Gamma];
    auto id = make_unit_matrix<std::complex<double>>(n);

    R_eq[Gamma][0] = -id;
    Rdot_eq[Gamma][0] = hamilt[Gamma][0];

    // first guess for R_eq
    for (int tau=1; tau<params.n_tau; tau++) R_eq[Gamma][tau] = 1.0;
  }

  // Initialize R_eq_old
  R_eq_old = R_eq;

  // solve Volterra equation
  std::vector<double> conv(n_blocks);
  conv[0] = 1;

  int nca_loops = 0;
  while (*std::max_element(conv.begin(), conv.end()) > params.tolerance) {

    nca_loops += 1;
    for (int Gamma=0; Gamma<n_blocks; Gamma++) {

      // get S_eq[Gamma] from R_eq
      for (int tau=0; tau<params.n_tau; tau++)
        get_S_from_NCA(Gamma, tau, 0, S_eq, R_eq, Delta_mat1, Delta_mat2);

      // solve Volterra equation
      int n = block_sizes[Gamma];
      matrix<std::complex<double>> H(n,n); H() = -1_j * hamilt[Gamma][0];
      for (int tau=1; tau<params.n_tau; tau++)
        volterra_step(R_eq[Gamma], Rdot_eq[Gamma], S_eq[Gamma], H, dtau, tau, 0);

    }

    for (int Gamma=0; Gamma<n_blocks; Gamma++) {
        conv[Gamma] = max_element(abs(R_eq[Gamma].data() - R_eq_old[Gamma].data())) / max_element(abs(R_eq[Gamma].data() + R_eq_old[Gamma].data()));
    }

    R_eq_old = R_eq;
  }

  // by convention S^mat = 1j*S^(0..-i\beta)
  for (int Gamma=0; Gamma<n_blocks; Gamma++)
    S_eq[Gamma] *= 1_j;

  std::cout << "Equilibrium part finished" << std::endl;
  std::cout << "Number of NCA loops: " << nca_loops << std::endl << std::endl;
}
