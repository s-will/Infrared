/**
 * @file libinfrared.cpp
 * @brief The libinfrared module exposing C++ infrared functionality to Python
 *
 * This file specifies the Python interface to the low level C++ library.
 * It is compiled into the Python module infrared.libinfrared. In Inrared's
 * high-level Python interface, most of the functionality is encapsulated
 * by Python classes and functions.
 */

/**
 * @package infrared.libinfrared
 * @copydoc libinfrared.cpp
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

namespace py = pybind11;
using namespace ired;

/**@brief trampoline / wrapper class
 *
 * @note this is needed to support overriding of virtual functions
 * in Function/Constraint objects from Python
 *
 * @private
 */
template<class FunValue>
struct PyFunction
    : public Function<FunValue> {
    using var_idx_t = int;
    using parent_t = Function<FunValue>;

    // inherit constructor
    //PyFunction(const std::vector<var_idx_t> &vars) : Function<FunValue>(vars) {}
    using parent_t::parent_t;

    FunValue
    operator () (const Assignment & a) const override {
        PYBIND11_OVERRIDE_PURE_NAME(
                FunValue,
                parent_t,
                "__call__",
                operator (),
                a);
    }

    std::string
    name() const override {
        PYBIND11_OVERRIDE(
            std::string,
            parent_t,
            name
            );
    }

    bool
    auto_materialize() const override {
        PYBIND11_OVERRIDE(
            bool,
            parent_t,
            auto_materialize
            );
    }
};

//! the type of partition functions
//using pf_type = long double;
using pf_type = __float128;

using PFClusterTree = ClusterTree<pf_type>;
using ArcticClusterTree = ClusterTree< int, ArcticEvaluationPolicy<int> >;


/**
 * @brief seed the C++-side random number generator
 */
void
seed(int x) {
    srand(x);
}


/** The libinfrared module exposing infrared to Python
 * @private
 */
PYBIND11_MODULE(libinfrared,ir)
{
    ir.doc() = "Infrared framework";
    ir.def("seed",&seed);

    py::class_< Assignment >
        (ir, "Assignment")
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

    py::class_< Function<pf_type>, PyFunction<pf_type>, std::shared_ptr<Function<pf_type>> >
        function(ir, "Function");

    function
        .def(py::init<const std::vector<int> &>())
        .def("__call__", &Function<pf_type>::operator ())
        .def("vars", &Function<pf_type>::vars,
             py::return_value_policy::copy)
        ;

    py::class_< Function<int>, PyFunction<int>, std::shared_ptr<Function<int>> >
        intfunction(ir, "IntFunction");

    intfunction
        .def(py::init<const std::vector<int> &>())
        .def("__call__", &Function<int>::operator ())
        .def("vars", &Function<int>::vars,
             py::return_value_policy::copy)
        ;

    py::class_< FiniteDomain >(ir, "FiniteDomain" )
        .def(py::init<int>())
        .def(py::init<int,int>())
        .def(py::init<std::pair<int,int>>())
        .def("lb", &FiniteDomain::lb)
        .def("ub", &FiniteDomain::ub)
        .def("empty", &FiniteDomain::empty)
        .def("size", &FiniteDomain::size)
        .def("undet", &FiniteDomain::undet)
        .def("in", &FiniteDomain::in)
        .def("contains", &FiniteDomain::in)
        .def("__copy__", [](const FiniteDomain &fd) {return FiniteDomain(fd);})
        .def("__deepcopy__",
                [](const FiniteDomain &fd, py::dict)
                  {return FiniteDomain(fd);}, "memo" )
        ;

    py::class_< PFClusterTree >(ir,"PFClusterTree")
        .def(py::init<int,FiniteDomain>())
        .def(py::init<const FiniteDomainVector &>())
        .def(py::init<int,int>())
        .def(py::init<const std::vector<int> &>())
        .def("add_root_cluster", &PFClusterTree::add_root_cluster)
        .def("add_child_cluster", &PFClusterTree::add_child_cluster)
        .def("add_constraint", &PFClusterTree::add_constraint)
        .def("add_function", &PFClusterTree::add_function)
        .def("evaluate", &PFClusterTree::evaluate)
        .def("is_consistent", &PFClusterTree::is_consistent)
        .def("resample", &PFClusterTree::restricted_traceback)
        .def("sample", &PFClusterTree::traceback)
        .def("sample_new", &PFClusterTree::traceback_new)
        .def("evaluation", &PFClusterTree::accumulated_weight)
        ;

    py::class_< ArcticClusterTree >(ir,"ArcticClusterTree")
        .def(py::init<int,FiniteDomain>())
        .def(py::init<const FiniteDomainVector &>())
        .def("add_root_cluster", &ArcticClusterTree::add_root_cluster)
        .def("add_child_cluster", &ArcticClusterTree::add_child_cluster)
        .def("add_constraint", &ArcticClusterTree::add_constraint)
        .def("add_function", &ArcticClusterTree::add_function)
        .def("evaluate", &ArcticClusterTree::evaluate)
        .def("is_consistent", &ArcticClusterTree::is_consistent)
        .def("optimize", &ArcticClusterTree::traceback)
        ;
}
