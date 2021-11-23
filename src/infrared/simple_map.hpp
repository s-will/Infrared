#ifndef SIMPLE_MAP_HH
#define SIMPLE_MAP_HH

/*
 * InfraRed ---  A generic engine for Boltzmann sampling over constraint networks
 * (C) Sebastian Will, 2018
 *
 * This file is part of the InfraRed source code.
 *
 * InfraRed provides a generic framework for tree decomposition-based
 * Boltzmann sampling over constraint networks
 */

/**
 * @file
 *
 * @brief Defines a simple (write-once) map class as space-saving replacement of hash_map
 */

#include <vector>
#include <algorithm>
#include <cassert>

namespace ired {

/**
 * @brief Space saving replacement for map
 *
 * Maintains lists of values in idx sorted lists, allowing
 * log-time access and efficient traversal. Erasing elements is
 * supported, but takes linear time.  Still, this data structure seems
 * to be a good compromise, since e.g. balanced trees or hashs require
 * a lot of space.
 */
template<class key_t, class val_t>
class SimpleMap: std::vector<std::pair<key_t, val_t> > {
    typedef std::pair<key_t, val_t> key_val_t;
    typedef std::vector<key_val_t> key_val_vec_t;
    typedef typename key_val_vec_t::iterator iterator;
    typedef typename key_val_vec_t::const_iterator const_iterator;

    static
    class {
    public:
	bool
	operator () (const key_val_t &x,
		     const key_t &y) const {
	    return x.first < y;
	}
    } comp;

    iterator
    binsearch (iterator first, iterator last, const key_t& key)
    {
	first = std::lower_bound(first,last,key,comp);
	if (first==last || key < first->first) {
	    return last;
	}
	return first;
    }

    const_iterator
    binsearch (const_iterator first, const_iterator last, const key_t& key) const
    {
	first = std::lower_bound(first,last,key,comp);
	if (first==last || key < first->first) {
	    return last;
	}
	return first;
    }

public:
    SimpleMap() {}

    const_iterator
    find(const key_t &key) const {
	auto it= binsearch(key_val_vec_t::begin(), key_val_vec_t::end(), key);
	assert(it == key_val_vec_t::end() || it->first == key);
	return it;
    };

    iterator
    find(const key_t &key) {
	auto it=binsearch(key_val_vec_t::begin(), key_val_vec_t::end(), key);
	assert(it == key_val_vec_t::end() || it->first == key);
	return it;
    };

    bool
    exists(const key_t &key) const {
	return find(key) != key_val_vec_t::end();
    }

    auto
    end() const {
        return key_val_vec_t::end();
    }

    /**
     * @brief push in ascending order of keys
     * @param key
     * @param val
     *
     * successive push_ascending must be in ascending order of the key type
     */
    void
    push_ascending( const key_t &key, const val_t &val ) {
	assert(size()==0||key > key_val_vec_t::operator[](size()-1).first);
	key_val_vec_t::push_back(key_val_t(key,val));
    }

    void
    erase(iterator it) {
	//std::cout << key_val_vec_t::size()<<" "<<key_val_vec_t::capacity()<<std::endl;

	// doing only this, wastes space (due to stl-vector allocation strategy)
	// key_val_vec_t::erase(it);

	// instead: copy to new space with exactly the right size
	// (possible optimization: only mark erased entries at
	// first and delay the copying until it pays off)
	key_val_vec_t vec(key_val_vec_t::size()-1);
	iterator target=copy(key_val_vec_t::begin(),it,vec.begin());
	copy(it+1,key_val_vec_t::end(),target);

	vec.swap(*this);
    }

    size_t
    size() const {
	return key_val_vec_t::size();
    }

    size_t
    capacity() const {
	return key_val_vec_t::capacity();
    }

    void
    reallocate() {
	key_val_vec_t vec(size());
	copy(key_val_vec_t::begin(),key_val_vec_t::end(),vec.begin());
	vec.swap(*this);
    }
};

} // end namespace ired

#endif //SIMPLE_MAP_HH
