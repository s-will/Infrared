#ifndef BOOST_PYTHON_AUX
#define BOOST_PYTHON_AUX

#include <boost/python.hpp>

namespace bpy = boost::python;

template<class T>
struct container_to_list
{
    static PyObject* convert(const T& xs)
    {
        bpy::list* l = new bpy::list();
        for(const auto &x: xs) {
            l->append(x);
        }

        return l->ptr();
    }
};

template<typename containedType>
struct convert_vector_from_pyseq{
	convert_vector_from_pyseq(){ 
            bpy::converter
                ::registry::push_back(&convertible,
                                      &construct,
                                      bpy::type_id<std::vector<containedType> >());
        }
	static void* convertible(PyObject* obj_ptr){
		// the second condition is important, for some reason otherwise there were attempted conversions of Body to list which failed afterwards.
		if(!PySequence_Check(obj_ptr) || !PyObject_HasAttrString(obj_ptr,"__len__")) return 0;
		return obj_ptr;
	}
	static void construct(PyObject* obj_ptr, bpy::converter::rvalue_from_python_stage1_data* data){
		 void* storage=((bpy::converter::rvalue_from_python_storage<std::vector<containedType> >*)(data))->storage.bytes;
		 new (storage) std::vector<containedType>();
		 std::vector<containedType>* v=(std::vector<containedType>*)(storage);
		 int l=PySequence_Size(obj_ptr); if(l<0) abort(); /*std::cerr<<"l="<<l<<"; "<<typeid(containedType).name()<<std::endl;*/ v->reserve(l); for(int i=0; i<l; i++) { v->push_back(bpy::extract<containedType>(PySequence_GetItem(obj_ptr,i))); }
		 data->convertible=storage;
	}
};


template<class T>
void
register_vector_conversions() {

    // c++ vector type that we want to convert to and from Python
    using register_t = std::vector<T>;

    const auto &info = bpy::type_id< register_t >(); 
    auto* reg = bpy::converter::registry::query(info); 
    if (reg == nullptr || (*reg).m_to_python == nullptr) {
        // install to Python converter
        bpy::to_python_converter< register_t , container_to_list<std::vector<T>> >();
    }

    // install from Python converter
    convert_vector_from_pyseq<T>();
}

#endif
