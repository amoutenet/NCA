#pragma once

#include <triqs/gfs.hpp>
#include <triqs/atom_diag/atom_diag.hpp>
#include "parameters.hpp"


using namespace triqs::gfs;

using gf_tau = gf<imtime>;
using block_gf_tau = block_gf<imtime>;
using block_gf_2tau = block_gf<cartesian_product<imtime, imtime>>;

using gf_t = gf<retime>;
using gf_2t = gf<cartesian_product<retime, retime>>;

using gf_ttau = gf<cartesian_product<retime, imtime>>;

template<typename M1, typename M2>
using block_gf_2t = block_gf<cartesian_product<M1, M2>>;

using mat = triqs::arrays::matrix<std::complex<double>>;


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

    // solver the equilibrium equations
    void solve_equilibrium();

    // solver the greater equations
    void solve_greater();

    // solver the mixed equations
    void solve_mixed();
    void get_mixed_Q(int Gamma, int t, int tau);

    // solver the lesser equations
    void solve_lesser();
    void get_lesser_Q(int Gamma, int t, int tp);
    void get_Rdot(int Gamma, int t, int tp);

    // get the Green functions
    void compute_G_mat();
    void compute_G_les();
    void compute_G_gtr();
    void compute_G_mix();

    // get S from R using the NCA equation - equilibrium overload
    void get_S_from_NCA(int Gamma, int t, int tp, block_gf_tau & S, block_gf_tau const & R, block_gf_2tau const & Delta_p, block_gf_2tau const & Delta_m);

    // get S from R using the NCA equation - mixed, lesser, greater overloads
    template<typename M1, typename M2>
    void get_S_from_NCA(int Gamma, int t, int tp, block_gf_2t<M1,M2> & S, block_gf_2t<M1,M2> const & R, block_gf_2t<M1,M2> const & Delta_p, block_gf_2t<M2,M1> const & Delta_m, int sign);

  
  public:

    triqs::hilbert_space::fundamental_operator_set fops;

    // parameters
    parameters params;
  
    // atomic diagonalization
    triqs::atom_diag::atom_diag<true> atom;
    // initialize atom diag
    void initialize_atom_diag(std::function<triqs::operators::many_body_operator_generic<std::complex<double>>(double)> function);
  
    // info for the blocks
    int n_blocks;
    std::vector<int> block_sizes;
    std::vector<int> n_particles;
  
    // discretization times
    double dt, dtau;
  
    // hybridization functions
    block_gf<cartesian_product<retime, retime>> Delta_gtr, Delta_les;
    block_gf<cartesian_product<retime, imtime>> Delta_mix;
    block_gf<cartesian_product<imtime, retime>> Delta_mix2;
    block_gf<imtime> Delta_mat;
    block_gf<cartesian_product<imtime, imtime>> Delta_mat1, Delta_mat2;
  
    // Green functions
    block_gf<imtime> G_mat;
    block_gf<cartesian_product<retime, retime>> G_gtr, G_les;
    block_gf<cartesian_product<retime, imtime>> G_mix;
  
    // hamiltonian matrix
    block_gf<retime> hamilt;
  
    // R and S matrices: S[\Gamma](Ntau, Ntau, n, n)
    block_gf<imtime> R_eq, Rdot_eq, S_eq;
    block_gf<cartesian_product<retime, retime>> R_gtr, Rdot_gtr, S_gtr,
                                                R_les, Rdot_les, S_les, Q_les;
    block_gf<cartesian_product<retime, imtime>> R_mix, Rdot_mix, S_mix, Q_mix;
  
  
    CPP2PY_ARG_AS_DICT
    solver(constructor_p const & cparams);
  
    // complete the special hybridization functions
    void complete_delta();
  
    // launch the solver
    CPP2PY_ARG_AS_DICT
    void solve(solve_p const & sparams);

    // get the partition function
    std::complex<double> get_Z();

    // make a single Volterra step
    // -- equilibrium overload
    void volterra_step(gf_tau & R, gf_tau & Rdot, gf_tau const & S, mat H, double dt, int t, int a);
    // -- greater and lesser overload
    void volterra_step(gf_2t & R, gf_2t & Rdot, gf_2t const & S, gf_t const & H, gf_2t const & Q, double dt, int t, int tp, int a);
    // -- mixed overload
    void volterra_step(gf_ttau & R, gf_ttau & Rdot, gf_2t const & S, gf_t const & H, gf_ttau const & Q, double dt, int t, int tp, int a);

};
