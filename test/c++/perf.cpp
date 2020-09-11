#include <iostream>
#include <triqs/gfs.hpp>
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>

#include <nca/solver.hpp>

using namespace triqs::gfs;
using namespace triqs::arrays;

using triqs::operators::c;
using triqs::operators::c_dag;
using triqs::operators::n;
using triqs::operators::many_body_operator;
using triqs::hilbert_space::gf_struct_t;

// ----------------------------------
// Kanamori Hamilt - performance check
// -----------------------------------

int main(int argc, char* argv[]) {
  triqs::mpi::communicator world;
  triqs::mpi::environment env(argc, argv);

 int beta = 2;
 double JoverU = 0.2;
 double U = 2;
 double J = JoverU * U;
 double Delta = 1.;
 double V = 1.;

 double mu = 2. * (U-2*J);

 double U_p = U - 2*J;
 double U_pp = U - 3*J;

 int nt = 2000;
 double tmax = 60;

 gf_struct_t gf_struct{{"up_0", {0}}, {"up_1", {0}}, {"up_2", {0}}, {"down_0", {0}}, {"down_1", {0}}, {"down_2", {0}}};

 // Construct solver
 solver nca_solver({gf_struct, tmax, nt});


 // Define and fill real-time Delta
 gf<retime> r_t{{-tmax,tmax, 2*nt-1}, {1,1}};
 gf<refreq> r_w = make_gf_from_fourier(r_t);
 for (auto &w : r_w.mesh()) {
   auto omega = real(dcomplex(w)); 
   if (std::abs(omega) < 10) r_w[w] = (0.1 * omega - 1_j * std::sqrt(1 - omega * omega / 100))/5;
   else r_w[w] = (0.1 * omega - (omega > 0 ? 1 : -1) * std::sqrt(omega * omega / 100 - 1))/5;
 }
 
 gf<refreq> Gamma_r = r_w;
 gf<refreq> Gamma_k = r_w;

 for (auto &w : Gamma_r.mesh()) {
   Gamma_r[w] = r_w[w];
 }
 
 for (auto const &w : Gamma_k.mesh()) {
   Gamma_k[w] = std::tanh(beta * real(dcomplex(w)) / 2.) * (Gamma_r[w] - conj(Gamma_r[w]));
 }
 
 gf<refreq> F_r = r_w;
 gf<refreq> F_k = r_w;

 for (auto &w : Gamma_r.mesh()) {
   F_r[w] = 2 * r_w[w];
 }
 
 for (auto const &w : F_k.mesh()) {
   F_k[w] = 0.5 * (std::tanh(beta * (real(dcomplex(w)) + V) / 2.) + std::tanh(beta * (real(dcomplex(w)) - V) / 2.)) * (F_r[w] - conj(F_r[w]));
 }


 gf<refreq> hyb_les_w = r_w;
 gf<refreq> hyb_gtr_w = r_w;

 for (auto const &w : hyb_les_w.mesh()) {
   auto fac = 2_j * imag(Gamma_r[w] + F_r[w]);
   hyb_les_w[w] = (Gamma_k[w] + F_k[w] - fac) / 2;
   hyb_gtr_w[w] = (Gamma_k[w] + F_k[w] + fac) / 2;
 }

 auto hyb_les_t = make_gf_from_fourier(hyb_les_w);
 auto hyb_gtr_t = make_gf_from_fourier(hyb_gtr_w);

 for (auto &delta: nca_solver.Delta_gtr) {
   for (auto const & t: delta.mesh()) {
     delta[t] = hyb_gtr_t(real(dcomplex(t)));
   }
 }
 
 for (auto &delta: nca_solver.Delta_les) {
   for (auto const & t: delta.mesh()) {
     delta[t] = hyb_les_t(real(dcomplex(t)));
   }
 }


 // Construct Kanamori Hamiltonian
 std::vector<double> chemical_pot {mu + Delta, mu, mu};
 std::vector<std::string> spin_names{"up","down"};

 many_body_operator H;

 // density terms
 double Uval = 0;
 for (auto const & s1 : spin_names) {
   for (auto const & s2 : spin_names) {
     for (int a1 = 0; a1 < 3; a1++) {
       for (int a2 =0; a2 < 3; a2++) {
	 if ((s1 == s2) and (a1 == a2)) Uval = -2 * chemical_pot[a1];
	 else if (s1 == s2) Uval = U_pp;
	 else if (a1 == a2) Uval = U;
	 else Uval = U_p;

	 H += 0.5 * Uval * n(s1,a1) * n(s2,a2);
       }
     }
   }
 }

 // spin-flip terms
 for (auto const & s1 : spin_names) {
   for (auto const & s2 : spin_names) {
     if (s1 == s2) continue;

     for (int a1 = 0; a1 < 3; a1++) {
       for (int a2 =0; a2 < 3; a2++) {
	 if (a1 == a2) continue;

	 H -= 0.5 * J * c_dag(s1,a1) * c(s2,a1) * c_dag(s2,a2) * c(s1,a2);
       }
     }
   }
 }

 // pair-hopping terms
 for (auto const & s1 : spin_names) {
   for (auto const & s2 : spin_names) {
     if (s1 == s2) continue;

     for (int a1 = 0; a1 < 3; a1++) {
       for (int a2 =0; a2 < 3; a2++) {
	 if (a1 == a2) continue;

	 H += 0.5 * J * c_dag(s1,a1) * c_dag(s2,a1) * c(s2,a2) * c(s1,a2);
       }
     }
   }
 }


 nca_solver.initialize_atom_diag(H);
 /*
   S.initialize_atom_diag(hamilt)

   R_init = S.R_gtr.copy()
   for Gamma in range(S.n_blocks):
          for it in range(nt-1, 2*nt-1):
	            for a in range(S.block_sizes[Gamma]):
		                  R_init['%i'%Gamma][MeshPoint(it)][a,a] = -1j
				              
				  p = {}
 p["H"] = hamilt
   p["R_init"] = R_init
#p["tolerance"] = 1e-5
#p["recompute_subspaces"] = False

   S.solve(**p)
 
   */

 return 0;
}




