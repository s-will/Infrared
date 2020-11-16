/**
 * @file
 * @brief The libinfrared module exposing infrared to Python
 *
 * This file specifies the python interface.
 */

/*
 * InfraRed ---  A generic engine for Boltzmann sampling over constraint networks
 * (C) Sebastian Will, 2018
 *
 * This file is part of the InfraRed source code.
 *
 * InfraRed provides a generic framework for tree decomposition-based
 * Boltzmann sampling over constraint networks
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <memory>

#include "infrared.hpp"
#include "rnadesign.hpp"

namespace py = pybind11; 
using namespace ired;
using namespace ired::rnadesign;

/**@brief trampoline / wrapper class
 *
 * @note this is needed to support overriding of the virtual function
 * '__call__' from Python
 */
template<class FunValue=double>
struct PyFunction
    : public Function<FunValue> {
    using var_idx_t = int;
    using parent_t = Function<FunValue>;  

    // inherit constructor
    //PyFunction(const std::vector<var_idx_t> &vars) : Function<FunValue>(vars) {}
    using parent_t::parent_t;

    FunValue
    operator () (const Assignment & a) const override { 
        PYBIND11_OVERLOAD_PURE_NAME(
                FunValue,
                parent_t,
                "__call__",
                operator (),
                a);
    }
};

/**
 * @brief seed the random number generator (for sampling, as well as
 * tree decomposition via libhtd)
 */
void
seed(int x) {
    srand(x);
}

//! @brief The libinfrared module exposing infrared to Python
PYBIND11_MODULE(libinfrared,ir)
{
    ir.doc() = "Infrared module or Boltzmann sampling";
    ir.def("seed",&seed);

    py::class_< Assignment >
        (ir, "Assignment")
        .def(py::init<int>())
        .def("values", &Assignment::values,
             py::return_value_policy::copy)
        ;

    py::class_< Constraint, PyFunction<bool>, std::shared_ptr<Constraint> >
        constraint(ir, "Constraint");
    
    constraint
        .def(py::init<const std::vector<int> &>())
        .def("__call__", &Constraint::operator ())
        .def("vars", &Constraint::vars,
             py::return_value_policy::copy)
        ;

    py::class_< Function<>, PyFunction<>, std::shared_ptr<Function<>> >
        function(ir, "Function");

    function
        .def(py::init<const std::vector<int> &>())
        //.def("__call__", &Function<>::operator ())
        .def("vars", &Function<>::vars,
             py::return_value_policy::copy)
        ;

    py::class_< ComplConstraint >
        (ir, "ComplConstraint", constraint)
        .def(py::init<int,int>())
       ;

    py::class_< SameComplClassConstraint >
        (ir, "SameComplClassConstraint", constraint)
        .def(py::init<int,int>())
       ;

    py::class_< DifferentComplClassConstraint >
        (ir, "DifferentComplClassConstraint", constraint )
        .def(py::init<int,int>())
       ;

    py::class_< BPEnergy >
        (ir, "BPEnergy", function)
        .def(py::init<int, int, bool, double>())
        .def("set_energy_table", &BPEnergy::set_energy_table)
        ;

    py::class_< StackEnergy >
        (ir, "StackEnergy", function)
        .def(py::init<int, int, double>())
        .def("set_energy_table", &StackEnergy::set_energy_table)
        ;

    py::class_< GCControl >
            (ir, "GCControl", function)
            .def(py::init<int, double>());

    py::class_< ClusterTree<> >(ir,"ClusterTree")
        .def(py::init<int,int>())
        .def(py::init<const std::vector<int> &>())
        .def("add_root_cluster", &ClusterTree<>::add_root_cluster)
        .def("add_child_cluster", &ClusterTree<>::add_child_cluster)
        .def("add_constraint", &ClusterTree<>::add_constraint)
        .def("add_function", &ClusterTree<>::add_function)
        .def("evaluate", &ClusterTree<>::evaluate)
        .def("is_consistent", &ClusterTree<>::is_consistent)
        .def("sample", &ClusterTree<>::sample)
        ;
}
