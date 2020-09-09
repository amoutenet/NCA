#pragma once

using namespace triqs::arrays;

// helper function
void solver::get_lesser_Q(int Gamma, int t, int tp) {

  Q_les[Gamma][{t,tp}] = 0.0;

  for (int k=0; k<tp; k++)
    Q_les[Gamma][{t,tp}] -= 0.5 * dt * (
         S_les[Gamma][{t,k}] * R_gtr[Gamma][{k,tp}] +
         S_les[Gamma][{t,k+1}] * R_gtr[Gamma][{k+1,tp}] );

  for (int k=0; k<params.n_tau-1; k++)
    Q_les[Gamma][{t,tp}] += 1_j * 0.5 * dtau * (
      S_mix[Gamma][{t,k}] * dagger(triqs::make_view(R_mix[Gamma][{tp,params.n_tau-1-k}])) +
      S_mix[Gamma][{t,k+1}] * dagger(triqs::make_view(R_mix[Gamma][{tp,params.n_tau-1-k-1}])) );

}


// helper function
void solver::get_Rdot(int Gamma, int t, int tp) {

  Rdot_les[Gamma][{t,tp}] = 0.0;

  for (int k=0; k<tp; k++)
    Rdot_les[Gamma][{t,tp}] -= 0.5 * dt * (
         S_les[Gamma][{t,k}] * R_gtr[Gamma][{k,tp}] +
         S_les[Gamma][{t,k+1}] * R_gtr[Gamma][{k+1,tp}] );

  for (int k=0; k<params.n_tau-1; k++)
    Rdot_les[Gamma][{t,tp}] += 1_j * 0.5 * dtau * (
      S_mix[Gamma][{t,k}] * dagger(triqs::make_view(R_mix[Gamma][{tp,params.n_tau-1-k}])) +
      S_mix[Gamma][{t,k+1}] * dagger(triqs::make_view(R_mix[Gamma][{tp,params.n_tau-1-k-1}])) );

  Rdot_les[Gamma][{t,tp}] += hamilt[Gamma][t] * R_les[Gamma][{t,tp}];

  for (int k=0; k<t; k++)
    Rdot_les[Gamma][{t,tp}] += 0.5 * dt * (
      S_gtr[Gamma][{t,k}] * R_les[Gamma][{k,tp}] +
      S_gtr[Gamma][{t,k+1}] * R_les[Gamma][{k+1,tp}] );

  Rdot_les[Gamma][{t,tp}] *= -1_j;

}



void solver::solve_lesser() {

  auto R_les_old = R_les;

  // first column t = 0
  for (int Gamma=0; Gamma<n_blocks; Gamma++)
    R_les[Gamma][{0,0}] = 1_j * R_eq[Gamma][params.n_tau-1];

  for (int Gamma=0; Gamma<n_blocks; Gamma++) {
    get_S_from_NCA(Gamma, 0, 0, S_les, R_les, Delta_les, Delta_gtr, -1);
    get_lesser_Q(Gamma, 0, 0);
    Rdot_les[Gamma][{0,0}] = -1_j * hamilt[Gamma][0] * R_les[Gamma][{0,0}] \
                               -1_j * Q_les[Gamma][{0,0}];
  }

  // Move in time
  int nca_loops = 0;
  for (int it=1; it<params.n_t; it++) {

    // Initial guess for R
    for (int Gamma=0; Gamma<n_blocks; Gamma++){
      for (int tp=0; tp<it+1; tp++) R_les[Gamma][{it,tp}] = 1.0;
    }
    R_les_old = R_les;

    // NCA loop
    std::vector<double> conv(n_blocks);
    conv[0] = 1;
    
    nca_loops = 0;
    while (*std::max_element(conv.begin(), conv.end()) > params.tolerance) {
      nca_loops += 1;

      // get R by solving Volterra equation for all blocks
      for (int Gamma=0; Gamma<n_blocks; Gamma++) {

        // get S_les[Gamma] from R_les
        for (int tp=0; tp<it; tp++) {
          get_S_from_NCA(Gamma, it, tp, S_les, R_les, Delta_les, Delta_gtr, -1);
          get_lesser_Q(Gamma, it, tp);
          volterra_step(R_les[Gamma], Rdot_les[Gamma], S_gtr[Gamma], hamilt[Gamma], Q_les[Gamma], dt, it, tp, 0);
        }
      }

      // Get R on the other side of the diagonal (all components needed for S below)
      for (int Gamma=0; Gamma<n_blocks; Gamma++) {
        for (int t=0; t<it; t++) R_les[Gamma][{t,it}] = -dagger(triqs::make_view(R_les[Gamma][{it,t}]));
      }

      // complete other side of diagonal
      for (int Gamma=0; Gamma<n_blocks; Gamma++) {

        for (int t=0; t<it; t++) {
          get_S_from_NCA(Gamma, t, it, S_les, R_les, Delta_les, Delta_gtr, -1);
          get_lesser_Q(Gamma, t, it);
          get_Rdot(Gamma, t, it);
        }

        // Diagonal element
        get_S_from_NCA(Gamma, it, it, S_les, R_les, Delta_les, Delta_gtr, -1);
        get_lesser_Q(Gamma, it, it);
        volterra_step(R_les[Gamma], Rdot_les[Gamma], S_gtr[Gamma], hamilt[Gamma], Q_les[Gamma], dt, it, it, 0);

      }

      for (int Gamma = 0; Gamma < n_blocks; Gamma++){
	conv[Gamma] = max_element(abs(R_les[Gamma].data() - R_les_old[Gamma].data())) / max_element(abs(R_les[Gamma].data() + R_les_old[Gamma].data()));
      }

      R_les_old = R_les;

    } // nca loop

  } // time loop

  std::cout << "Lesser part finished" << std::endl;
  std::cout << "Number of NCA loops for the last time point: " << nca_loops << std::endl << std::endl;

}
