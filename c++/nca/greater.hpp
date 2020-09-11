#pragma once

#include <numeric>

using namespace triqs::arrays;

void solver::solve_greater() {

  // initial condition and guess
  for (int Gamma=0; Gamma<n_blocks; Gamma++) {
    int n = block_sizes[Gamma];
    auto id = make_unit_matrix<std::complex<double>>(n);

    // first guess for R_gtr
    R_gtr[Gamma] = params.R_init[Gamma];
    for (int it = 0; it<params.n_t - 1; it++) S_tilde_gtr[Gamma][it] = 0 * id;
  }

  // Initialize R_gtr_old
  auto R_gtr_old = R_gtr;

  // NCA loop
  std::vector<double> conv_abs(n_blocks);
  std::vector<double> conv_rel(n_blocks);
  double err_rel = 1.;

  int nca_loops = 0;
  while ((err_rel > params.tolerance) && (nca_loops < 200)) {
    nca_loops +=1;

    for (int Gamma=0; Gamma<n_blocks; Gamma++) {
      int n = block_sizes[Gamma];
      auto id = make_unit_matrix<std::complex<double>>(n);
      auto km = make_zero_tail(R_gtr_w[Gamma], 4);

      // get S_gtr[Gamma] from R_gtr
      get_S_from_NCA(Gamma, S_tilde_gtr, R_gtr, Delta_gtr, Delta_les, 1);
      //for (int it = 0; it<params.n_t-1; it++) S_tilde_gtr[Gamma][it] = 0 * id;

      make_ft_S_tilde(Gamma);

      // get R_gtr[Gamma] from S_gtr
      for (auto const & w : R_gtr_w[Gamma].mesh()) {
	auto x = imag(inverse(real(dcomplex(w))* id - hamilt[Gamma] - S_tilde_gtr_w[Gamma][w]));
	R_gtr_w[Gamma][w] = 2_j * x;
      }

      R_gtr[Gamma] = make_gf_from_fourier(R_gtr_w[Gamma], R_gtr[Gamma].mesh(), km);
      
      // Normalisation
      R_gtr[Gamma] = R_gtr[Gamma] * (-1_j * n) / trace(R_gtr[Gamma][params.n_t-1]);
    }

    for (int Gamma=0; Gamma < n_blocks; Gamma++) { 
      conv_abs[Gamma] = max_element(abs(R_gtr[Gamma].data() - R_gtr_old[Gamma].data()));
      
      if (conv_abs[Gamma] < 1e-9) conv_rel[Gamma] = 0.;
      else conv_rel[Gamma] = max_element(abs(R_gtr[Gamma].data() - R_gtr_old[Gamma].data())) / max_element(abs(R_gtr[Gamma].data() + R_gtr_old[Gamma].data()));
    }

    err_rel = *std::max_element(conv_rel.begin(), conv_rel.end());
    R_gtr_old = R_gtr;

  } // nca loop

  // Complete t < 0
  S_gtr = S_tilde_gtr;
  for (int Gamma=0; Gamma < n_blocks; Gamma++) {
    for (int it = 0; it<params.n_t - 1; it++){
      S_gtr[Gamma][it] = -dagger(triqs::make_view(S_gtr[Gamma][2*params.n_t-2-it]));
    }
  }

  // Get S_gtr_w
  for (int Gamma=0; Gamma<n_blocks; Gamma++) {
    for (auto const & w : S_gtr_w[Gamma].mesh()) {
      auto x = imag(S_tilde_gtr_w[Gamma][w]);
      S_gtr_w[Gamma][w] = 2_j * x;
    }
  }

  R_tilde_gtr = R_gtr;

  for (int Gamma=0; Gamma < n_blocks; Gamma++) {
    int n = block_sizes[Gamma];
    auto id = make_unit_matrix<std::complex<double>>(n);

    for (int it = params.n_t; it<2*params.n_t-1; it++) R_tilde_gtr[Gamma][it] = 0 * id;
  }

  make_ft_R_tilde();

  std::cout << "Greater part finished" << std::endl;
  std::cout << "Number of NCA loops: " << nca_loops << std::endl << std::endl;

}

void solver::make_ft_R_tilde() {

  for (int Gamma=0; Gamma < n_blocks; Gamma++) {

    for (int it = params.n_t; it < 2*params.n_t-1; it++) fit_R[Gamma][it]() = 0;

    int n = block_sizes[Gamma];
    for (int a=0; a<n; a++) {
      for (int b=0; b<n; b++) {

        for (int it = 0; it < params.n_t; it++) {
          rho[it] = abs(R_tilde_gtr[Gamma][it](a,b));
        }

        double norm = std::accumulate(rho.begin(), rho.end(), 0.0);

        double var = 0;
        for (int it = 0; it < params.n_t; it++) {
          double x = (it - params.n_t + 1) * dt;
          var += x * x * rho[it] / norm;
        }
        double stdev = std::sqrt(var);
        double Gam = 4./stdev;

        double r_0 = imag(R_tilde_gtr[Gamma][params.n_t - 1](a,b));
        double rdot_0 = -real(R_tilde_gtr[Gamma][params.n_t - 2](a,b)) / dt;

	if (abs(r_0) < 1e-9) {
          for (int it = 0; it < params.n_t; it++) {
	    fit_R[Gamma][it](a,b) = 0.;
	  }

          for (auto const & w : fit_R_w[Gamma].mesh()) {
	    fit_R_w[Gamma][w](a,b) = 0.;
	  }
	}

	else {
          double frac = rdot_0/r_0;

          double eps1 = (abs(frac) > 0.5) ? frac : frac + 2;
          double eps2 = (abs(frac) > 0.5) ? frac : frac + 1;

          for (int it = 0; it < params.n_t; it++) {
            double t = (it-params.n_t+1)*dt;
            fit_R[Gamma][it](a,b) = (-r_0) * 1_j * (std::exp(Gam*t - 1_j*eps1*t) - 2 * std::exp(0.5*Gam*t - 1_j*eps2*t));
          }

          for (auto const & w : fit_R_w[Gamma].mesh()) {
            fit_R_w[Gamma][w](a,b) = (-r_0) / (real(dcomplex(w)) - eps1 - 1_j*Gam) + 2*r_0 / (real(dcomplex(w)) - eps2 - 0.5*1_j*Gam);
          }
	}

      }
    }

    temp_R[Gamma] = R_tilde_gtr[Gamma] - fit_R[Gamma];

    auto km = make_zero_tail(temp_R_w[Gamma], 4);
    temp_R_w[Gamma] = make_gf_from_fourier(temp_R[Gamma], temp_R_w[Gamma].mesh(), km);

    R_tilde_gtr_w[Gamma] = temp_R_w[Gamma] + fit_R_w[Gamma];
    
  }

}


void solver::make_ft_S_tilde(int Gamma) {


    for (int it = 0; it < params.n_t-1; it++) fit_S[Gamma][it]() = 0;

    int n = block_sizes[Gamma];
    for (int a=0; a<n; a++) {
      for (int b=0; b<n; b++) {

        for (int it = 0; it < params.n_t; it++) {
          rho[it] = abs(S_tilde_gtr[Gamma][it+params.n_t-1](a,b));
        }

        double norm = std::accumulate(rho.begin(), rho.end(), 0.0);

        double var = 0;
        for (int it = 0; it < params.n_t; it++) {
          double x = it * dt;
          var += x * x * rho[it] / norm;
        }
        double stdev = std::sqrt(var);
        double Gam = 4./stdev;

        double s_0 = imag(S_tilde_gtr[Gamma][params.n_t - 1](a,b));
        double sdot_0 = real(S_tilde_gtr[Gamma][params.n_t](a,b)) / dt;

	if (abs(s_0) < 1e-9) {
          for (int it = params.n_t - 1; it < 2*params.n_t-1; it++) {
            fit_S[Gamma][it](a,b) = 0.;
          }

          for (auto const & w : fit_S_w[Gamma].mesh()) {
            fit_S_w[Gamma][w](a,b) = 0.;
          }
	}

	else {
          double frac = sdot_0/s_0;

          double eps1 = (abs(frac) > 0.5) ? frac : frac + 2;
          double eps2 = (abs(frac) > 0.5) ? frac : frac + 1;


          for (int it = params.n_t - 1; it < 2*params.n_t-1; it++) {
            double t = (it-params.n_t+1)*dt;
            fit_S[Gamma][it](a,b) = (-s_0) * 1_j * (std::exp(-Gam*t - 1_j*eps1*t) - 2 * std::exp(-0.5*Gam*t - 1_j*eps2*t));
          }

          for (auto const & w : fit_S_w[Gamma].mesh()) {
            fit_S_w[Gamma][w](a,b) = s_0 / (real(dcomplex(w)) - eps1 + 1_j*Gam) - 2*s_0 / (real(dcomplex(w)) - eps2 + 0.5 *1_j*Gam);
          }
	}

      }
    }

    temp_S[Gamma] = S_tilde_gtr[Gamma] - fit_S[Gamma];

    auto km = make_zero_tail(temp_S_w[Gamma], 4);
    temp_S_w[Gamma] = make_gf_from_fourier(temp_S[Gamma], temp_S_w[Gamma].mesh(), km);

    S_tilde_gtr_w[Gamma] = temp_S_w[Gamma] + fit_S_w[Gamma];

}
