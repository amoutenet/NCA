#pragma once

#include <triqs/arrays/blas_lapack/gemm.hpp>
#include <triqs/arrays/linalg/det_and_inverse.hpp>

/*

Below is the Volterra solver that, given the matrix-valued functions $S(t,t')$
and $H(t)$ returns the function $R(t,t')$ that satisfies

$$ i \frac{\partial R(t,t')}{\partial t} = H(t) R(t,t')
   + \int_{a}^t S(t,t'') R(t'',t') \mathrm{d}t'' + Q(t,t') $$

with different initial conditions.

*/


using namespace triqs::arrays;

using gf_tau = triqs::gfs::gf<triqs::gfs::imtime>;
using mat = matrix<std::complex<double>>;

// Equilibrium version overload -- Q is always 0
void solver::volterra_step(gf_tau & R, gf_tau & Rdot, gf_tau const & S, mat H, double dt, int t, int a) {

  int n = R.target_shape()[0];

  // identity
  auto id = make_unit_matrix<std::complex<double>>(n);
  matrix<std::complex<double>> M(n,n);
  M() = 0;

  // convolution S*R without last time
  for (int j=a+1; j<t; j++) M += S[t-j] * R[j];
  M += 0.5 * S[t-a] * R[a];

  // compute R
  auto A = inverse(id + 1_j * (dt/2.) * (H + (dt/2.)*S[0]));
  auto B = R[t-1] + (dt/2.) * (Rdot[t-1] - 1_j * dt * M);
  R[t] = A * B;

  // complete convolution with last time
  M += 0.5 * S[0] * R[t];

  // compute Rdot
  Rdot[t] = -1_j * H * R[t] - 1_j * dt * M;

}


using gf_t = triqs::gfs::gf<triqs::gfs::retime>;
using gf_2t = triqs::gfs::gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::retime>>;

// Greater and lesser version overload
void solver::volterra_step(gf_2t & R, gf_2t & Rdot, gf_2t const & S, gf_t const & H, gf_2t const & Q, double dt, int t, int tp, int a) {

  int n = R.target_shape()[0];

  // identity
  auto id = make_unit_matrix<std::complex<double>>(n);
  matrix<std::complex<double>> M(n,n);
  M() = 0;

  // convolution S*R without last time
  for (int j=a+1; j<t; j++) M += S[{t,j}] * R[{j,tp}];
  M += 0.5 * S[{t,a}] * R[{a,tp}];

  // compute R
  auto A = inverse(id + 1_j * (dt/2.) * (H[t] + (dt/2.)*S[{t,t}]));
  auto B = R[{t-1,tp}] + (dt/2.) * (Rdot[{t-1,tp}] - 1_j * Q[{t,tp}] - 1_j * dt * M);
  R[{t,tp}] = A * B;

  // complete convolution with last time
  M += 0.5 * S[{t,t}] * R[{t,tp}];

  // compute Rdot
  Rdot[{t,tp}] = -1_j * H[t] * R[{t,tp}] - 1_j * Q[{t,tp}] - 1_j * dt * M;

}


using gf_ttau = triqs::gfs::gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::imtime>>;

// Mixed version overload
void solver::volterra_step(gf_ttau & R, gf_ttau & Rdot, gf_2t const & S, gf_t const & H, gf_ttau const & Q, double dt, int t, int tp, int a) {

  int n = R.target_shape()[0];

  // identity
  auto id = make_unit_matrix<std::complex<double>>(n);
  matrix<std::complex<double>> M(n,n);
  M() = 0;

  // convolution S*R without last time
  for (int j=a+1; j<t; j++) M += S[{t,j}] * R[{j,tp}];
  M += 0.5 * S[{t,a}] * R[{a,tp}];

  // compute R
  auto A = inverse(id + 1_j * (dt/2.) * (H[t] + (dt/2.)*S[{t,t}]));
  auto B = R[{t-1,tp}] + (dt/2.) * (Rdot[{t-1,tp}] - 1_j * Q[{t,tp}] - 1_j * dt * M);
  R[{t,tp}] = A * B;

  // complete convolution with last time
  M += 0.5 * S[{t,t}] * R[{t,tp}];

  // compute Rdot
  Rdot[{t,tp}] = -1_j * H[t] * R[{t,tp}] - 1_j * Q[{t,tp}] - 1_j * dt * M;

}
