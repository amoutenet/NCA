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

using aview = triqs::arrays::array_view<std::complex<double>, 3>;
using mview = triqs::arrays::matrix_view<std::complex<double>>;

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

c.add_member(c_name = "dt",
             c_type = "double",
             read_only= False,
             doc = """""")

c.add_member(c_name = "dtau",
             c_type = "double",
             read_only= False,
             doc = """""")

c.add_member(c_name = "Delta_gtr",
             c_type = "triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::retime> >",
             read_only= False,
             doc = """""")

c.add_member(c_name = "Delta_les",
             c_type = "triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::retime> >",
             read_only= False,
             doc = """""")

c.add_member(c_name = "Delta_mix",
             c_type = "triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::imtime> >",
             read_only= False,
             doc = """""")

c.add_member(c_name = "Delta_mix2",
             c_type = "triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::imtime, triqs::gfs::retime> >",
             read_only= False,
             doc = """""")

c.add_member(c_name = "Delta_mat",
             c_type = "triqs::gfs::block_gf<triqs::gfs::imtime>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "Delta_mat1",
             c_type = "triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::imtime, triqs::gfs::imtime> >",
             read_only= False,
             doc = """""")

c.add_member(c_name = "Delta_mat2",
             c_type = "triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::imtime, triqs::gfs::imtime> >",
             read_only= False,
             doc = """""")

c.add_member(c_name = "G_mat",
             c_type = "triqs::gfs::block_gf<triqs::gfs::imtime>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "G_gtr",
             c_type = "triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::retime> >",
             read_only= False,
             doc = """""")

c.add_member(c_name = "G_les",
             c_type = "triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::retime> >",
             read_only= False,
             doc = """""")

c.add_member(c_name = "G_mix",
             c_type = "triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::imtime> >",
             read_only= False,
             doc = """""")

c.add_member(c_name = "hamilt",
             c_type = "triqs::gfs::block_gf<triqs::gfs::retime>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "R_eq",
             c_type = "triqs::gfs::block_gf<triqs::gfs::imtime>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "Rdot_eq",
             c_type = "triqs::gfs::block_gf<triqs::gfs::imtime>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "S_eq",
             c_type = "triqs::gfs::block_gf<triqs::gfs::imtime>",
             read_only= False,
             doc = """""")

c.add_member(c_name = "R_gtr",
             c_type = "triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::retime> >",
             read_only= False,
             doc = """""")

c.add_member(c_name = "Rdot_gtr",
             c_type = "triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::retime> >",
             read_only= False,
             doc = """""")

c.add_member(c_name = "S_gtr",
             c_type = "triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::retime> >",
             read_only= False,
             doc = """""")

c.add_member(c_name = "R_les",
             c_type = "triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::retime> >",
             read_only= False,
             doc = """""")

c.add_member(c_name = "Rdot_les",
             c_type = "triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::retime> >",
             read_only= False,
             doc = """""")

c.add_member(c_name = "S_les",
             c_type = "triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::retime> >",
             read_only= False,
             doc = """""")

c.add_member(c_name = "Q_les",
             c_type = "triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::retime> >",
             read_only= False,
             doc = """""")

c.add_member(c_name = "R_mix",
             c_type = "triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::imtime> >",
             read_only= False,
             doc = """""")

c.add_member(c_name = "Rdot_mix",
             c_type = "triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::imtime> >",
             read_only= False,
             doc = """""")

c.add_member(c_name = "S_mix",
             c_type = "triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::imtime> >",
             read_only= False,
             doc = """""")

c.add_member(c_name = "Q_mix",
             c_type = "triqs::gfs::block_gf<triqs::gfs::cartesian_product<triqs::gfs::retime, triqs::gfs::imtime> >",
             read_only= False,
             doc = """""")

c.add_constructor("""(**constructor_p)""", doc = """
+----------------+-------------------------------------------------------------------------------------------+---------+---------------+
| Parameter Name | Type                                                                                      | Default | Documentation |
+================+===========================================================================================+=========+===============+
| gf_struct      | triqs::hilbert_space::gf_struct_t                                                         |         |               |
+----------------+-------------------------------------------------------------------------------------------+---------+---------------+
| beta           | double                                                                                    |         |               |
+----------------+-------------------------------------------------------------------------------------------+---------+---------------+
| n_tau          | int                                                                                       |         |               |
+----------------+-------------------------------------------------------------------------------------------+---------+---------------+
| t_max          | double                                                                                    |         |               |
+----------------+-------------------------------------------------------------------------------------------+---------+---------------+
| n_t            | int                                                                                       |         |               |
+----------------+-------------------------------------------------------------------------------------------+---------+---------------+
""")

c.add_method("""void solve (**solve_p)""",
             doc = """
+---------------------+-------------------------------------------------------------------------------------------+---------+---------------+
| Parameter Name      | Type                                                                                      | Default | Documentation |
+=====================+===========================================================================================+=========+===============+
| hamilt              | std::function<triqs::operators::many_body_operator_generic<std::complex<double>>(double)> |         |               |
+---------------------+-------------------------------------------------------------------------------------------+---------+---------------+
| tolerance           | double                                                                                    | 1e-6    |               |
+---------------------+-------------------------------------------------------------------------------------------+---------+---------------+
| eq_solver           | bool                                                                                      | false   |               |
+---------------------+-------------------------------------------------------------------------------------------+---------+---------------+
| recompute_subspaces | bool                                                                                      | true    |               |
+---------------------+-------------------------------------------------------------------------------------------+---------+---------------+
""")

c.add_method("std::complex<double> get_Z()",
               doc = """""")

module.add_class(c)


# Converter for solve_p
c = converter_(
        c_type = "solve_p",
        doc = """""",
)

c.add_member(c_name = "hamilt",
    c_type = "std::function<triqs::operators::many_body_operator_generic<std::complex<double>>(double)>",
             initializer = """  """,
             doc = """""")

c.add_member(c_name = "tolerance",
             c_type = "double",
             initializer = """ 1e-6 """,
             doc = """""")

c.add_member(c_name = "eq_solver",
             c_type = "bool",
             initializer = """ false """,
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

c.add_member(c_name = "beta",
             c_type = "double",
             initializer = """  """,
             doc = """""")

c.add_member(c_name = "n_tau",
             c_type = "int",
             initializer = """ """,
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
