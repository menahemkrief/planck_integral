#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "planck_integral.hpp"

namespace planck_integral {
    void bind_planck_integral(pybind11::module& m){
        using namespace pybind11::literals;
        m.def("planck_integral", &planck_integral, pybind11::kw_only(), "a"_a, "b"_a);
        m.def("planck_energy_density_group_integral", &planck_energy_density_group_integral, pybind11::kw_only(), "E_low"_a, "E_high"_a, "T"_a);
    }
}

PYBIND11_MODULE(_planck_integral, m){
    m.doc() = "planck integral c++ module";

    planck_integral::bind_planck_integral(m);
}