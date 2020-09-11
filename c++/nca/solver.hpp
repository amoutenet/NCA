#pragma once

#include <triqs/gfs.hpp>
#include <triqs/atom_diag/atom_diag.hpp>
#include "parameters.hpp"

using gf_t = triqs::gfs::gf<triqs::gfs::retime>;
using block_gf_t = triqs::gfs::block_gf<triqs::gfs::retime>;

using gf_2t = triqs::gfs::gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::retime>>;
using block_gf_2t = triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::retime>>;


struct solver {

  private:

    // communicator
    triqs::mpi::communicator world;

    // update hamiltonian
    void update_hamiltonian(std::function<triqs::operators::many_body_operator_generic<std::complex<double>>(double)> function);

    // initialize all the different hybridization functions
    void initialize_hybridization_functions();
    // initialize all R and S matrices
    void initialize_R_and_S();
    // initialize hamiltonian
    void initialize_hamiltonian();

    // solver the greater equations
    void solve_greater();
    // solver the lesser equations
    void solve_lesser(std::vector<triqs::arrays::array<std::complex<double>,2>> const & R_init);
    void get_lesser_Q(int Gamma, int t, int tp);
    void get_Rdot(int Gamma, int t, int tp);

    // get the Green functions
    void compute_G_les();
    void compute_G_gtr();

    // Compute S from R
    void get_S_from_NCA(int Gamma, int t, int tp, block_gf_2t & S, block_gf_2t const & R, block_gf_2t const & Delta_p, block_gf_2t const & Delta_m, int sign);


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
    block_gf_2t Delta_gtr, Delta_les;

    // Green functions
    block_gf_2t G_gtr, G_les;

    // hamiltonian matrix
    block_gf_t hamilt;

    // R and S matrices: S[\Gamma](Ntau, Ntau, n, n)
    block_gf_2t R_gtr, Rdot_gtr, S_gtr,
                R_les, Rdot_les, S_les, Q_les;

    CPP2PY_ARG_AS_DICT
    solver(constructor_p const & cparams);

    // initialize atom diag
    void initialize_atom_diag(std::function<triqs::operators::many_body_operator_generic<std::complex<double>>(double)> function);

    // launch the solver
    CPP2PY_ARG_AS_DICT
    void solve(solve_p const & sparams);

    // get the partition function
    std::complex<double> get_Z();

    // make a single Volterra step
    void volterra_step(gf_2t & R, gf_2t & Rdot, gf_2t const & S, gf_t const & H, gf_2t const & Q, double dt, int t, int tp, int a);
};
