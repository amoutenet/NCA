#pragma once

using namespace triqs::arrays;
using mview = matrix_view<std::complex<double>>;

std::complex<double> solver::get_Z() {

  std::complex<double> Z = 0;

  for (int Gamma=0; Gamma<n_blocks; Gamma++)
      Z += 1_j * trace(R_les[Gamma][params.n_t-1]);
  return Z;

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

            for (int it=0; it<2*params.n_t-1; it++)
                G_les[bln][it](a, b) -= 1_j * trace(R_gtr[Gamma][2*params.n_t-2-it] *
                                            c_mat * R_les[Gamma_p][it] * cdag_mat);

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

            for (int it=0; it<2*params.n_t-1; it++)
                G_gtr[bln][it](a, b) += 1_j * trace(R_les[Gamma][2*params.n_t-2-it] *
                                            c_mat * R_gtr[Gamma_p][it] * cdag_mat);

         }
        }
      }
    }
    bln++; // increase gf block number
  }

  G_gtr /= get_Z();

}
