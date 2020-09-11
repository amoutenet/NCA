#pragma once

using namespace triqs::arrays;
using block_gf_t = triqs::gfs::block_gf<triqs::gfs::retime>;

void solver::get_S_from_NCA(int Gamma, block_gf_t & S, block_gf_t const & R, block_gf_t const & Delta_p, block_gf_t const & Delta_m, int sign) {

  for (int it = params.n_t-1; it < 2*params.n_t-1; it++) {
    S[Gamma][it] = 0;

    // First term
    int bln = 0;
    for (auto const & [bname, idx_list] : params.gf_struct) {
      auto bl_size = idx_list.size();
      for (int a=0; a<bl_size; a++) {
        for (int b=0; b<bl_size; b++) {

          int ind_c = atom.get_fops()[{bname, b}];
          int ind_cdag = atom.get_fops()[{bname, a}];

          int Gamma_p = atom.c_connection(ind_c, Gamma);
          if (Gamma_p != -1) {

            TRIQS_ASSERT(atom.cdag_connection(ind_cdag, Gamma_p) == Gamma);
            auto c_mat = atom.c_matrix(ind_c, Gamma);
            auto cdag_mat = atom.cdag_matrix(ind_cdag, Gamma_p);

            S[Gamma][it] += 1_j * sign * cdag_mat * R[Gamma_p][it] * c_mat * Delta_p[bln][it](a,b);

          }
        }
      }
      bln++; // increase gf block number
    }

    // Second term
    bln = 0;
    for (auto const & [bname, idx_list] : params.gf_struct) {
      auto bl_size = idx_list.size();
      for (int a=0; a<bl_size; a++) {
        for (int b=0; b<bl_size; b++) {

          int ind_c = atom.get_fops()[{bname, b}];
          int ind_cdag = atom.get_fops()[{bname, a}];

          int Gamma_p = atom.cdag_connection(ind_cdag, Gamma);
          if (Gamma_p != -1) {

            TRIQS_ASSERT(atom.c_connection(ind_c, Gamma_p) == Gamma);
            auto c_mat = atom.c_matrix(ind_c, Gamma_p);
            auto cdag_mat = atom.cdag_matrix(ind_cdag, Gamma);

            S[Gamma][it] -= 1_j * sign * c_mat * R[Gamma_p][it] * cdag_mat * Delta_m[bln][2*params.n_t-2-it](a,b);

          }
        }
      }
      bln++; // increase gf block number
    }

  }

}
