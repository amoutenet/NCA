#pragma once

using namespace triqs::arrays;

void solver::solve_greater() {

  auto _ = range();
  auto R_gtr_old = R_gtr;
  auto Q = R_gtr; 

  for (int Gamma=0; Gamma<n_blocks; Gamma++) {
    int n = block_sizes[Gamma];
    auto id = make_unit_matrix<std::complex<double>>(n);

    for (auto const & [t,tp] : Q[Gamma].mesh()) {
       Q[Gamma][t,tp] = 0. * id;
    }
  } 


  // initial condition and guess
  for (int Gamma=0; Gamma<n_blocks; Gamma++) {

    int n = block_sizes[Gamma];
    auto id = make_unit_matrix<std::complex<double>>(n);

    R_gtr[Gamma][{0,0}] = -1_j * id;
    Rdot_gtr[Gamma][{0,0}] = -hamilt[Gamma][0];

  }

  // get S at (0,0) from R
  for (int Gamma=0; Gamma<n_blocks; Gamma++)
    get_S_from_NCA(Gamma, 0, 0, S_gtr, R_gtr, Delta_gtr, Delta_les, 1);


  // t > 0
  int nca_loops = 0;
  for (int it=1; it<params.n_t; it++) {

    // first guess for R_gtr
    for (int Gamma=0; Gamma<n_blocks; Gamma++) {
      int n = block_sizes[Gamma];
      auto id = make_unit_matrix<std::complex<double>>(n);

      for (int tp=0; tp<it+1; tp++) R_gtr[Gamma][{it,tp}] = 1.0 * id;
    }

    // Initialize R_gtr_old
    R_gtr_old = R_gtr;

    // NCA loop
    std::vector<double> conv(n_blocks);
    conv[0] = 1;
    
    nca_loops = 0;
    while (*std::max_element(conv.begin(), conv.end()) > params.tolerance) {
      nca_loops +=1;

      // get R by solving Volterra equation for all blocks
      for (int Gamma=0; Gamma<n_blocks; Gamma++) {

        // get S_gtr[Gamma] from R_gtr
        for (int tp=0; tp<it+1; tp++) {
          get_S_from_NCA(Gamma, it, tp, S_gtr, R_gtr, Delta_gtr, Delta_les,1);
          if (tp < it) get_S_from_NCA(Gamma, tp, it, S_gtr, R_gtr, Delta_gtr, Delta_les,1);
        }

	int n = block_sizes[Gamma];

        for (int j=0; j<it; j++) {
          volterra_step(R_gtr[Gamma], Rdot_gtr[Gamma], S_gtr[Gamma], hamilt[Gamma], Q[Gamma], dt, it, j, j);
          R_gtr[Gamma][{j,it}] = -dagger(triqs::make_view(R_gtr[Gamma][{it,j}]));
        }

        // diagonal element
        auto id = make_unit_matrix<std::complex<double>>(n);
        R_gtr[Gamma][{it,it}] = -1_j * id;
        Rdot_gtr[Gamma][{it,it}] = -hamilt[Gamma][it];

      }

      for (int Gamma=0; Gamma < n_blocks; Gamma++) {
	conv[Gamma] = max_element(abs(R_gtr[Gamma].data() - R_gtr_old[Gamma].data())) / max_element(abs(R_gtr[Gamma].data() + R_gtr_old[Gamma].data()));
      }

      R_gtr_old = R_gtr;

    } // nca loop

  } // t>0 loop

  std::cout << "Greater part finished" << std::endl;
  std::cout << "Number of NCA loops for last time point: " << nca_loops << std::endl << std::endl;

}
