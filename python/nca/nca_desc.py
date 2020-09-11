# Generated automatically using the command :
# c++2py ../../c++/nca/solver.hpp -p --cxxflags="-std=c++17" -C pytriqs --only="solver parameters solve_p constructor_p"
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "nca", doc = "", app_name = "nca")

# Imports
import pytriqs.atom_diag
import pytriqs.gf
import pytriqs.operators

# Add here all includes
module.add_include("nca/solver.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/complex.hpp>
#include <cpp2py/converters/pair.hpp>
#include <cpp2py/converters/string.hpp>
#include <cpp2py/converters/variant.hpp>
#include <cpp2py/converters/vector.hpp>
#include <triqs/cpp2py_converters/arrays.hpp>
#include <triqs/cpp2py_converters/gf.hpp>
#include <triqs/cpp2py_converters/operators_real_complex.hpp>

""")


# The class solver
c = class_(
        py_type = "Solver",  # name of the python class
        c_type = "solver",   # name of the C++ class
        doc = """""",   # doc of the C++ class
        hdf5 = False,
)

c.add_member(c_name = "atom",
             c_type = "triqs::atom_diag::atom_diag<true>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "n_blocks",
             c_type = "int",
             read_only= False,
             doc = """""")

c.add_member(c_name = "block_sizes",
             c_type = "std::vector<int>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "n_particles",
             c_type = "std::vector<int>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "parity",
             c_type = "std::vector<int>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "dt",
             c_type = "double",
             read_only= False,
             doc = """""")

c.add_member(c_name = "Delta_gtr",
             c_type = "triqs::gfs::block_gf<triqs::gfs::retime>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "Delta_les",
             c_type = "triqs::gfs::block_gf<triqs::gfs::retime>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "G_gtr",
             c_type = "triqs::gfs::block_gf<triqs::gfs::retime>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "G_les",
             c_type = "triqs::gfs::block_gf<triqs::gfs::retime>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "hamilt",
        c_type = "std::vector<triqs::arrays::matrix<std::complex<double>>>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "R_gtr",
             c_type = "triqs::gfs::block_gf<triqs::gfs::retime>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "R_tilde_gtr",
             c_type = "triqs::gfs::block_gf<triqs::gfs::retime>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "S_gtr",
             c_type = "triqs::gfs::block_gf<triqs::gfs::retime>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "S_tilde_gtr",
             c_type = "triqs::gfs::block_gf<triqs::gfs::retime>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "R_les",
             c_type = "triqs::gfs::block_gf<triqs::gfs::retime>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "S_les",
             c_type = "triqs::gfs::block_gf<triqs::gfs::retime>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "temp_R",
             c_type = "triqs::gfs::block_gf<triqs::gfs::retime>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "fit_R",
             c_type = "triqs::gfs::block_gf<triqs::gfs::retime>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "temp_S",
             c_type = "triqs::gfs::block_gf<triqs::gfs::retime>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "fit_S",
             c_type = "triqs::gfs::block_gf<triqs::gfs::retime>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "R_gtr_w",
             c_type = "triqs::gfs::block_gf<triqs::gfs::refreq>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "R_tilde_gtr_w",
             c_type = "triqs::gfs::block_gf<triqs::gfs::refreq>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "S_tilde_gtr_w",
             c_type = "triqs::gfs::block_gf<triqs::gfs::refreq>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "R_les_w",
             c_type = "triqs::gfs::block_gf<triqs::gfs::refreq>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "S_les_w",
             c_type = "triqs::gfs::block_gf<triqs::gfs::refreq>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "temp_R_w",
             c_type = "triqs::gfs::block_gf<triqs::gfs::refreq>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "fit_R_w",
             c_type = "triqs::gfs::block_gf<triqs::gfs::refreq>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "temp_S_w",
             c_type = "triqs::gfs::block_gf<triqs::gfs::refreq>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "fit_S_w",
             c_type = "triqs::gfs::block_gf<triqs::gfs::refreq>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "rho",
        c_type = "std::vector<double>",
             read_only= False,
             doc = """""")


c.add_constructor("""(**constructor_p)""", doc = """
+----------------+-------------------------------------------------------------------------------------------+---------+---------------+
| Parameter Name | Type                                                                                      | Default | Documentation |
+================+===========================================================================================+=========+===============+
| gf_struct      | triqs::hilbert_space::gf_struct_t                                                         |         |               |
+----------------+-------------------------------------------------------------------------------------------+---------+---------------+
| t_max          | double                                                                                    |         |               |
+----------------+-------------------------------------------------------------------------------------------+---------+---------------+
| n_t            | int                                                                                       |         |               |
+----------------+-------------------------------------------------------------------------------------------+---------+---------------+
""")

c.add_method("""void solve (**solve_p)""",
             doc = """
+--------------------------+-------------------------------------------------------------------------------------------+---------+---------------+
| Parameter Name           | Type                                                                                      | Default | Documentation |
+==========================+===========================================================================================+=========+===============+
| H                        | triqs::operators::many_body_operator_generic<std::complex<double>>                        |         |               |
+--------------------------+-------------------------------------------------------------------------------------------+---------+---------------+
| R_init                   | triqs::gfs::block_gf<triqs::gfs::retime>                                                  |         |               |
+--------------------------+-------------------------------------------------------------------------------------------+---------+---------------+
| tolerance                | double                                                                                    | 1e-6    |               |
+--------------------------+-------------------------------------------------------------------------------------------+---------+---------------+
| S_threshold              | double                                                                                    | 1e-9    |               |
+--------------------------+-------------------------------------------------------------------------------------------+---------+---------------+
| recompute_subspaces      | bool                                                                                      | true    |               |
+--------------------------+-------------------------------------------------------------------------------------------+---------+---------------+
""")

c.add_method("void initialize_atom_diag (triqs::operators::many_body_operator_generic<std::complex<double>> H)",
               doc = """""")

c.add_method("std::complex<double> get_Z()",
               doc = """""")

module.add_class(c)


# Converter for solve_p
c = converter_(
        c_type = "solve_p",
        doc = """""",
)

c.add_member(c_name = "H",
    c_type = "triqs::operators::many_body_operator_generic<std::complex<double>>",
             initializer = """  """,
             doc = """""")

c.add_member(c_name = "R_init",
             c_type = "triqs::gfs::block_gf<triqs::gfs::retime>",
             initializer = """  """,
             doc = """""")

c.add_member(c_name = "tolerance",
             c_type = "double",
             initializer = """ 1e-6 """,
             doc = """""")

c.add_member(c_name = "S_threshold",
             c_type = "double",
             initializer = """ 1e-9 """,
             doc = """""")

c.add_member(c_name = "recompute_subspaces",
             c_type = "bool",
             initializer = """ true """,
             doc = """""")

module.add_converter(c)

# Converter for constructor_p
c = converter_(
        c_type = "constructor_p",
        doc = """""",
)
c.add_member(c_name = "gf_struct",
             c_type = "triqs::hilbert_space::gf_struct_t",
             initializer = """  """,
             doc = """""")

c.add_member(c_name = "t_max",
             c_type = "double",
             initializer = """  """,
             doc = """""")

c.add_member(c_name = "n_t",
             c_type = "int",
             initializer = """  """,
             doc = """""")


module.add_converter(c)


module.generate_code()
