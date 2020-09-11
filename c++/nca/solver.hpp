#pragma once

#include <triqs/gfs.hpp>
#include <triqs/atom_diag/atom_diag.hpp>
#include "parameters.hpp"

using gf_t = triqs::gfs::gf<triqs::gfs::retime>;
using block_gf_t = triqs::gfs::block_gf<triqs::gfs::retime>;
using block_gf_w = triqs::gfs::block_gf<triqs::gfs::refreq>;
using mat = triqs::arrays::matrix<std::complex<double>>;

struct solver {

  private:

    // communicator
    triqs::mpi::communicator world;

    // update hamiltonian
    void update_hamiltonian(triqs::operators::many_body_operator_generic<std::complex<double>> H);

    // initialize all the different hybridization functions
    void initialize_hybridization_functions();
    // initialize all R and S matrices
    void initialize_R_and_S();
    // initialize hamiltonian
    void initialize_hamiltonian();

    // solver the greater equations
    void solve_greater();
    // solver the lesser equations
    void solve_lesser();

    // get the Green functions
    void compute_G_les();
    void compute_G_gtr();

    // Compute S from R
    void get_S_from_NCA(int Gamma, block_gf_t & S, block_gf_t const & R, block_gf_t const & Delta_p, block_gf_t const & Delta_m, int sign);

    // Make Fourier transforms of R tilde and S tilde
    void make_ft_R_tilde();
    void make_ft_S_tilde(int Gamma);


  public:

    triqs::hilbert_space::fundamental_operator_set fops;

    // parameters
    parameters params;

    // atomic diagonalization
    triqs::atom_diag::atom_diag<true> atom;

    // info for the blocks
    int n_blocks;
    std::vector<int> block_sizes;
    std::vector<int> n_particles;
    std::vector<int> parity;

    // discretization times
    double dt;

    // hybridization functions
    block_gf_t Delta_gtr, Delta_les;

    // Green functions
    block_gf_t G_gtr, G_les;

    // hamiltonian matrix
    std::vector<mat> hamilt;

    // R and S matrices: S[\Gamma](Ntau, Ntau, n, n)
    block_gf_t R_gtr, S_gtr, R_tilde_gtr, S_tilde_gtr,
               R_les, S_les;
    block_gf_w R_tilde_gtr_w, S_tilde_gtr_w,
	       R_gtr_w, S_gtr_w, R_les_w, S_les_w;

    // Helpers for the FT
    block_gf_t fit_R, temp_R, fit_S, temp_S;
    block_gf_w fit_R_w, temp_R_w, fit_S_w, temp_S_w;
    std::vector<double> rho;

    CPP2PY_ARG_AS_DICT
    solver(constructor_p const & cparams);

    // initialize atom diag
    void initialize_atom_diag(triqs::operators::many_body_operator_generic<std::complex<double>> H);

    // launch the solver
    CPP2PY_ARG_AS_DICT
    void solve(solve_p const & sparams);

    // get the partition function
    std::complex<double> get_Z();
};
