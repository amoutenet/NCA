#pragma once

using namespace triqs::arrays;
using mview = matrix_view<std::complex<double>>;

std::complex<double> solver::get_Z() {

  std::complex<double> Z = 0;

  for (int Gamma=0; Gamma<n_blocks; Gamma++)
      Z -= trace(R_eq[Gamma][params.n_tau-1]);
  return Z;

}


void solver::compute_G_mat() {

  G_mat() = 0;

  int bln = 0;
  for (auto const & [bname, idx_lst] : params.gf_struct) {
    auto bl_size = idx_lst.size();
    for (int a=0; a<bl_size; a++) {
      for (int b=0; b<bl_size; b++) {

        int ind_c = atom.get_fops()[{bname, a}];
        int ind_cdag = atom.get_fops()[{bname, b}];

        for (int Gamma=0; Gamma<n_blocks; Gamma++) {
          int Gamma_p = atom.cdag_connection(ind_cdag, Gamma);
          if (Gamma_p != -1) {

            TRIQS_ASSERT(atom.c_connection(ind_c, Gamma_p) == Gamma);
            auto c_mat = atom.c_matrix(ind_c, Gamma_p);
            auto cdag_mat = atom.cdag_matrix(ind_cdag, Gamma);

            for (int itau=0; itau<params.n_tau; itau++)
              G_mat[bln][itau](a, b) -= trace(R_eq[Gamma][params.n_tau-1-itau] *
                                              c_mat * R_eq[Gamma_p][itau] * cdag_mat);

         }
        }
      }
    }
    bln++; // increase gf block number
  }

  G_mat /= get_Z();

}

void solver::compute_G_les() {

  G_les() = 0;

  int bln = 0;
  for (auto const & [bname, idx_lst] : params.gf_struct) {
    auto bl_size = idx_lst.size();
    for (int a=0; a<bl_size; a++) {
      for (int b=0; b<bl_size; b++) {

        int ind_c = atom.get_fops()[{bname, a}];
        int ind_cdag = atom.get_fops()[{bname, b}];

        for (int Gamma=0; Gamma<n_blocks; Gamma++) {
          int Gamma_p = atom.cdag_connection(ind_cdag, Gamma);
          if (Gamma_p != -1) {

            TRIQS_ASSERT(atom.c_connection(ind_c, Gamma_p) == Gamma);
            auto c_mat = atom.c_matrix(ind_c, Gamma_p);
            auto cdag_mat = atom.cdag_matrix(ind_cdag, Gamma);

            for (int t=0; t<params.n_t; t++)
              for (int tp=0; tp<params.n_t; tp++)
                G_les[bln][{t,tp}](a, b) -= 1_j * trace(R_gtr[Gamma][{tp,t}] *
                                            c_mat * R_les[Gamma_p][{t,tp}] * cdag_mat);

         }
        }
      }
    }
    bln++; // increase gf block number
  }

  G_les /= get_Z();

}

void solver::compute_G_gtr() {

  G_gtr() = 0;

  int bln = 0;
  for (auto const & [bname, idx_lst] : params.gf_struct) {
    auto bl_size = idx_lst.size();
    for (int a=0; a<bl_size; a++) {
      for (int b=0; b<bl_size; b++) {

        int ind_c = atom.get_fops()[{bname, a}];
        int ind_cdag = atom.get_fops()[{bname, b}];

        for (int Gamma=0; Gamma<n_blocks; Gamma++) {
          int Gamma_p = atom.cdag_connection(ind_cdag, Gamma);
          if (Gamma_p != -1) {

            TRIQS_ASSERT(atom.c_connection(ind_c, Gamma_p) == Gamma);
            auto c_mat = atom.c_matrix(ind_c, Gamma_p);
            auto cdag_mat = atom.cdag_matrix(ind_cdag, Gamma);

            for (int t=0; t<params.n_t; t++)
              for (int tp=0; tp<params.n_t; tp++)
                G_gtr[bln][{t,tp}](a, b) += 1_j * trace(R_les[Gamma][{tp,t}] *
                                            c_mat * R_gtr[Gamma_p][{t,tp}] * cdag_mat);

         }
        }
      }
    }
    bln++; // increase gf block number
  }

  G_gtr /= get_Z();

}

void solver::compute_G_mix() {

  G_mix() = 0;

  int bln = 0;
  for (auto const & [bname, idx_lst] : params.gf_struct) {
    auto bl_size = idx_lst.size();
    for (int a=0; a<bl_size; a++) {
      for (int b=0; b<bl_size; b++) {

        int ind_c = atom.get_fops()[{bname, a}];
        int ind_cdag = atom.get_fops()[{bname, b}];

        for (int Gamma=0; Gamma<n_blocks; Gamma++) {
          int Gamma_p = atom.cdag_connection(ind_cdag, Gamma);
          if (Gamma_p != -1) {

            TRIQS_ASSERT(atom.c_connection(ind_c, Gamma_p) == Gamma);
            auto c_mat = atom.c_matrix(ind_c, Gamma_p);
            auto cdag_mat = atom.cdag_matrix(ind_cdag, Gamma);

            for (int t=0; t<params.n_t; t++)
              for (int tau=0; tau<params.n_tau; tau++)
                G_mix[bln][{t,tau}](a, b) += 1_j * trace(dagger(triqs::make_view(R_mix[Gamma][{t,params.n_tau-1-tau}])) *
                                            c_mat * R_mix[Gamma_p][{t,tau}] * cdag_mat);

         }
        }
      }
    }
    bln++; // increase gf block number
  }

  G_mix /= get_Z();

}
