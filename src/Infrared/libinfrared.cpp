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

#include <boost/python.hpp>

#include "infrared.hpp"
#include "rnadesign.hpp"

#include <vector>

#include "boost_python_aux.hpp"

namespace bpy = boost::python;

using namespace ired;
using namespace ired::rnadesign;

template<class FunValue=double>
struct FunctionWrap
    : public Function<FunValue>, bpy::wrapper<Function<FunValue>> {
    using var_idx_t = int;

    FunctionWrap(const std::vector<var_idx_t> &vars)
        : Function<FunValue>(vars) {}

    FunValue
    operator () (const Assignment & a) const {
        return this->get_override("__call__")(a);
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
BOOST_PYTHON_MODULE(libinfrared)
{
    register_vector_conversions<int>();
    register_vector_conversions<double>();

    bpy::def("seed",&seed);

    bpy::class_<Assignment>("Assignment", bpy::init<int>())
        .def("values", &Assignment::values,
             bpy::return_value_policy<bpy::copy_const_reference>())
        ;

    bpy::class_< FunctionWrap<bool>,
            std::shared_ptr<FunctionWrap<bool>>,
            boost::noncopyable >
        ("Constraint", bpy::init<const std::vector<int> &>())
        .def("__call__", bpy::pure_virtual(&Constraint::operator ()))
        .def("vars", &Constraint::vars,
             bpy::return_value_policy<bpy::copy_const_reference>())
        ;

    bpy::class_< FunctionWrap<>,
            std::shared_ptr<FunctionWrap<>>,
            boost::noncopyable >
        ("Function", bpy::init<const std::vector<int> &>())
        .def("__call__", bpy::pure_virtual(&Function<>::operator ()))
        .def("vars", &Function<>::vars,
             bpy::return_value_policy<bpy::copy_const_reference>())
        ;

    bpy::implicitly_convertible< std::shared_ptr<FunctionWrap<>>,
                                 std::shared_ptr<Function<>>      >();

    bpy::class_<ComplConstraint, bpy::bases<Constraint>,
           boost::noncopyable>("ComplConstraint", bpy::init<int,int>())
       ;

    bpy::class_<SameComplClassConstraint, bpy::bases<Constraint>,
           boost::noncopyable>("SameComplClassConstraint", bpy::init<int,int>())
       ;

    bpy::class_<DifferentComplClassConstraint, bpy::bases<Constraint>,
           boost::noncopyable>("DifferentComplClassConstraint", bpy::init<int,int>())
       ;

    bpy::class_<BPEnergy, std::shared_ptr<BPEnergy>,
                bpy::bases<Function<>>, boost::noncopyable>
        ("BPEnergy", bpy::init<int, int, bool, double>())
        .def("set_energy_table", &BPEnergy::set_energy_table)
        ;

    bpy::class_<StackEnergy, std::shared_ptr<StackEnergy>,
                bpy::bases<Function<>>, boost::noncopyable>
        ("StackEnergy", bpy::init<int, int, double>())
        .def("set_energy_table", &StackEnergy::set_energy_table)
        ;

    bpy::class_<GCControl, bpy::bases<Function<>>,
           boost::noncopyable>("GCControl", bpy::init<int, double>());

    bpy::class_<ClusterTree<>>("ClusterTree", bpy::init<int,int>())
        .def("add_root_cluster", &ClusterTree<>::add_root_cluster)
        .def("add_child_cluster", &ClusterTree<>::add_child_cluster)
        .def("add_constraint", &ClusterTree<>::add_constraint)
        .def("add_function", &ClusterTree<>::add_function)
        .def("evaluate", &ClusterTree<>::evaluate)
        .def("sample", &ClusterTree<>::sample)
        ;
}
