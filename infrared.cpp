/*
 * InfraRed ---  A generic engine for Boltzmann sampling over constraint networks
 * (C) Sebastian Will, 2018
 *
 * This file is part of the InfraRed source code.
 *
 * InfraRed provides a generic framework for tree decomposition-based
 * Boltzmann sampling over constraint networks
 *
 * This file specifies the python interface.
 */


#include <boost/python.hpp>

#include "infrared.hpp"
#include "rnadesign.hpp"

#include <memory>

using namespace boost::python;

using namespace ired;

using boost::python::object;

template<class T>
struct vector_to_list
{
    static PyObject* convert(const std::vector<T>& vec)
    {
        boost::python::list* l = new boost::python::list();
        for(size_t i = 0; i < vec.size(); i++) {
            l->append(vec[i]);
        }

        return l->ptr();
    }
};


template<typename containedType>
struct convert_vector_from_seq{
	convert_vector_from_seq(){ converter::registry::push_back(&convertible,&construct,type_id<std::vector<containedType> >()); }
	static void* convertible(PyObject* obj_ptr){
		// the second condition is important, for some reason otherwise there were attempted conversions of Body to list which failed afterwards.
		if(!PySequence_Check(obj_ptr) || !PyObject_HasAttrString(obj_ptr,"__len__")) return 0;
		return obj_ptr;
	}
	static void construct(PyObject* obj_ptr, converter::rvalue_from_python_stage1_data* data){
		 void* storage=((converter::rvalue_from_python_storage<std::vector<containedType> >*)(data))->storage.bytes;
		 new (storage) std::vector<containedType>();
		 std::vector<containedType>* v=(std::vector<containedType>*)(storage);
		 int l=PySequence_Size(obj_ptr); if(l<0) abort(); /*std::cerr<<"l="<<l<<"; "<<typeid(containedType).name()<<std::endl;*/ v->reserve(l); for(int i=0; i<l; i++) { v->push_back(extract<containedType>(PySequence_GetItem(obj_ptr,i))); }
		 data->convertible=storage;
	}
};

template<class FunValue=double>
struct FunctionWrap
    : public Function<FunValue>, wrapper<Function<FunValue>> {
    using var_idx_t = int;

    FunctionWrap(const std::vector<var_idx_t> &vars)
        : Function<FunValue>(vars) {}

    FunValue
    operator () (const Assignment & a) const {
        return this->get_override("__call__")(a);
    }
};

void
seed(int x) {
    srand(x);
}



BOOST_PYTHON_MODULE(infrared)
{
    to_python_converter<std::vector<int>, vector_to_list<int> >();
    convert_vector_from_seq<int>();

    to_python_converter<std::vector<double>, vector_to_list<double> >();
    convert_vector_from_seq<double>();

    def("seed",&seed);

    class_<Assignment>("Assignment", init<int>())
        .def("values", &Assignment::values,
             return_value_policy<copy_const_reference>())
        .def_readonly("undetermined", &Assignment::Undetermined)
        ;

    class_< FunctionWrap<bool>,
            std::shared_ptr<FunctionWrap<bool>>,
            boost::noncopyable >
        ("Constraint", init<const std::vector<int> &>())
        .def("__call__", pure_virtual(&Constraint::operator ()))
        ;

    class_< FunctionWrap<>,
            std::shared_ptr<FunctionWrap<>>,
            boost::noncopyable >
        ("Function", init<const std::vector<int> &>())
        .def("__call__", pure_virtual(&Function<>::operator ()))
        ;

    implicitly_convertible< std::shared_ptr<FunctionWrap<>>,
                            std::shared_ptr<Function<>>      >();

    class_<ComplConstraint, bases<Constraint>,
           boost::noncopyable>("ComplConstraint", init<int,int>())
       ;

    class_<BPEnergy, std::shared_ptr<BPEnergy>, bases<Function<>>, boost::noncopyable>
        ("BPEnergy", init<int,int, double>())
        .def("set_energy_table", &BPEnergy::set_energy_table)
        ;

    class_<GCControl, bases<Function<>>,
           boost::noncopyable>("GCControl", init<int, double>());

    class_<ClusterTree<>>("ClusterTree", init<int,int>())
        .def("add_root_cluster", &ClusterTree<>::add_root_cluster)
        .def("add_child_cluster", &ClusterTree<>::add_child_cluster)
        .def("add_constraint", &ClusterTree<>::add_constraint)
        .def("add_function", &ClusterTree<>::add_function)
        .def("evaluate", &ClusterTree<>::evaluate)
        .def("sample", &ClusterTree<>::sample)
        ;

    // implicitly_convertible< std::shared_ptr<FunctionWrap<>>,
    //                         std::shared_ptr<Function<>>      >();

    // implicitly_convertible< std::shared_ptr<FunctionWrap<bool>>,
    //                         std::shared_ptr<Function<bool>>      >();

}
