#ifndef INFRARED_FUNCTIONS_HPP
#define INFRARED_FUNCTIONS_HPP

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
 * @brief Defines function and constraint base classes
 *
 * Custom functions and constraints inherit from the abstract classes
 * Functions<FunVal> and Constraint and specialize by overriding (at
 * least) operator ().
 *
 * MaterializedFunction offers functions that tabulate their
 * values. This is used for messages in the cluster tree
 * evaluation and pre-computing of function values.
 */

#include <vector>
#include <unordered_map>
#include <memory>
#include <string>
#include <iterator>
#include <cassert>
#include <numeric>

#include "assignment.hpp"
#include "simple_map.hpp"

namespace ired {

    template<class FunValue, class EP>
    class FeatureNetwork;

    /**
     * @brief Dependencies specify a dependency between variables
     */
    class Dependency {
    public:
        using var_idx_t = int;

        explicit
        Dependency(const std::vector<var_idx_t> &vars) : vars_(vars) {}

        const std::vector<var_idx_t> &
        vars() const {return vars_;}

        virtual
        ~Dependency() {}

    private:
        const std::vector<var_idx_t> vars_;
    };

    /**
     * @brief Functions evaluate assignments of a subset of variables
     */
    template<class FunValue=double>
    class Function : public Dependency {
    public:
        using self_t = Function<FunValue>;
        using base_t = self_t;
        using assignment_t = Assignment;
        using fun_value_t = FunValue;

        explicit
        Function(const std::vector<var_idx_t> &vars) : Dependency(vars) {}

        virtual
        fun_value_t
        operator () (const assignment_t &) const = 0;

        /** quick check whether function is definitely zero at assignment
         * @param a assignment
         * @return true, if certainly zero
         *
         * some functions know how to quick check for zero (in particular sparse materialized functions)
         */
        virtual
        bool
        guaranteed_zero(const assignment_t &a) const { return false; }

        virtual
        bool
        auto_materialize() const { return true; }

        virtual
        std::string
        name() const { return "Function"; }

        virtual
        ~Function() {}
    };

    /**
     * @brief WeightedFunction operates in the same way as @see Function
     * but with added constant weight to the value
     */
    template <class FunValue = double>
    class WeightedFunction : public Function<FunValue> {
    public:
        using base_t = Function<FunValue>;
        using function_ptr_t = const base_t*;
        using assignment_t = typename base_t::assignment_t;
        using fun_value_t = typename base_t::fun_value_t;
        using var_idx_t = int;
        using vertex_descriptor_t = size_t;

        WeightedFunction(const std::vector<var_idx_t> &vars,
                         function_ptr_t function, 
                         std::function<fun_value_t(vertex_descriptor_t, const assignment_t&)> lambda_function, 
                         vertex_descriptor_t vertex
        ) : base_t(vars),
            function(function),
            lambda_function(lambda_function),
            vertex(vertex) {}

        fun_value_t
        operator() (const assignment_t& a) const override {
            auto c = lambda_function(vertex, a);
            return (*function)(a) - c;
        }
    
    private:
        function_ptr_t function;
        std::function<fun_value_t(vertex_descriptor_t, const assignment_t&)> lambda_function;
        vertex_descriptor_t vertex;
    };

    template <class T>
    inline
    std::ostream & operator << (std::ostream &out, const Function<T> &f) {
        out << f.name() << "(";
        for (auto i : f.vars()) { out << i << ", "; }
        out << ") <"<<&f<<">";
        return out;
    }

    /**
     * @brief map selector class
     *
     * for selecting a sparse table to represent
     * materialized functions
     */
    struct mapS {};

    /**
     * @brief simple map selector class
     *
     * for selecting an extra space saving, specialized sparse table to represent
     * materialized functions
     */
    struct simple_mapS {};

    /**
     * @brief vector selector class
     *
     * for selecting a non-sparse table to represent
     * materialized functions
     */
    struct vecS {};

    /**
     * @brief Switching containers in MaterializedFuntion
     */
    template<class FunValue, class T>
    struct container_selector {
        using type = T;
    };

    /**
     * @brief Class implementing the specializations to enable
     * instantiation of non-sparse materialized function classes
     */
    template<class FunValue>
    struct container_selector<FunValue,vecS> {
        //using type = std::vector<FunValue>;
        using type = typename vector_nbv_sel<FunValue>::type;
        static void init(type &x, size_t size, const FunValue &zero) {
            x.resize(size);
            std::fill(x.begin(),x.end(),zero);
        }
        static bool guaranteed_zero(const type &x, size_t i) {
            return false;
        }
        static FunValue get(const type &x, size_t i, const FunValue &zero) {
            assert(i<x.size());
            return x[i];
        }
        static void set(type &x,size_t i, const FunValue &v, const FunValue &zero) {x[i]=v;}
    };

    /**
     * @brief Class implementing the specializations to enable
     * instantiation of sparse materialized function classes
     */
    template<class FunValue>
    struct container_selector<FunValue,mapS> {
        using type = std::unordered_map<int,FunValue>;
        static void init(type &x, size_t size, const FunValue &zero) {
        }
        static const FunValue& get(const type &x, size_t i, const FunValue &zero) {
            auto it = x.find(i);
            if (it != x.end()) return it->second; else return zero;
        }
        static bool guaranteed_zero(const type &x, size_t i) {
            return x.find(i) == x.end();
        }
        static void set(type &x,size_t i, const FunValue &v, const FunValue &zero) {
            if (v!=zero) {
                x[i]=v;
            } else {
                x.erase(i); //erases key i if it exists
            }
        }
    };

    /**
     * @brief Class implementing the specializations to enable
     * space saving instantiation of sparse materialized function classes
     * using the specialized 'map' implementation simple_map
     * @note this requires strictly in order write access and has log
     * read access time
     */
    template<class FunValue>
    struct container_selector<FunValue,simple_mapS> {
        using type = SimpleMap<int,FunValue>;
        static void init(type &x, size_t size, const FunValue &zero) {
        }
        static const FunValue& get(const type &x, size_t i, const FunValue &zero) {
            auto it = x.find(i);
            if (it != x.end()) return it->second; else return zero;
        }
        static bool guaranteed_zero(const type &x, size_t i) {
            return x.find(i) == x.end();
        }
        static void set( type &x, size_t i, const FunValue &v, const FunValue &zero ) {
            if (v!=zero) {
                x.push_ascending(i,v);
            } else {
                // x.erase(i); //erases key i if it exists
            }
        }
    };



    /**
     * @brief A materialized function
     *
     * A function that holds a table of function values. It is used to
     * represent computed messages (DP matrices) during evaluation of
     * the cluster tree, or precomputed functions.
     *
     * Supports sparse and non-sparse representation of the function
     * value table. The choice is made by a template mechanism using
     * selector classes.
     */
    template< class FunValue, class ContainerS = simple_mapS >
    class MaterializedFunction : public Function<FunValue> {
    public:
        using self_t = MaterializedFunction<FunValue,ContainerS>;
        using parent_t = Function<FunValue>;
        using base_t = typename parent_t::base_t;

        using var_idx_t = typename parent_t::var_idx_t;
        using assignment_t = typename parent_t::assignment_t;
        using fun_value_t = FunValue;
        using function_t = Function<fun_value_t>;

        using data_t = typename container_selector<FunValue,ContainerS>::type;

        /**
         * @brief construct 'empty' with variables, domains, and zero value
         *
         * @param vars vector of indices of variable in a cn
         * @param domains domains of the variables in the cn
         * @param zero value (default for not explicitly set function values)
         *
         * The (non-zero) function values are typically set after
         * construction, using method @see set()
         */
        template<class FeatureNetwork>
        MaterializedFunction(const std::vector<var_idx_t> &vars,
                             const FeatureNetwork &cn
                             )
            :
            parent_t(vars),
            domains_(extract_domains(cn.domains())),
            zero_(FeatureNetwork::evaluation_policy_t::zero()),
            name_("Message")
        {
            container_selector<FunValue,ContainerS>::init(data_, calc_size(), zero_);
        }

        /**
         * @brief materializing constructor
         *
         * copy constructs from function, thereby materializing it
         *
         * @param function the function to be copied
         * @param cn constraint network (must outlive this object)
         * @param evaluation policy
         *
         */
        template<class FeatureNetwork>
        MaterializedFunction(const function_t *function,
                             const FeatureNetwork &cn
                             )
            :
            parent_t(unique(function->vars())),
            domains_(extract_domains(cn.domains())),
            zero_(FeatureNetwork::evaluation_policy_t::zero()),
            name_(function->name())
        {
            using feature_network_t = FeatureNetwork;

            using ep = typename feature_network_t::evaluation_policy_t;

            container_selector<FunValue,ContainerS>::init(data_, calc_size(), zero_);

            auto a = Assignment(cn.domains());

            auto constraints = std::vector<const typename feature_network_t::constraint_t *>();
            auto functions = std::vector<const typename feature_network_t::function_t *> {function};

            auto initial_value = a.eval_determined(functions, ep()); //evaluate 0-ary functions

            auto it = a.make_iterator(this->vars(),
                                      cn,
                                      constraints,
                                      functions,
                                      initial_value
                                      );

            for(; ! it.finished() ; ++it ) {
                set( a, it.value() );
            }
        }

        /**
         * @brief evaluate function
         *
         * @param a assignment where the function is evaluated
         */
        fun_value_t
        operator () ( const assignment_t & a ) const override {
            return container_selector<FunValue,ContainerS>::get(data_, index(a), zero_);
        }

        /**
         * @brief set function value
         *
         * @param a assignment where the value is set
         * @param value value to be set
         */
        void
        set( const assignment_t & a, const fun_value_t &val) {
            return container_selector<FunValue,ContainerS>::set(data_, index(a), val, zero_);
        }

        /**
         * @brief check whether this object knows that some function value is zero
         *
         * @param a assignment where this is checked
         */
        bool guaranteed_zero(const assignment_t &a) const override {
            return container_selector<FunValue,ContainerS>::guaranteed_zero(data_, index(a));
        }

        //! @brief whether to automatically materialize when added to CN
        virtual
        bool
        auto_materialize() const override { return false; }

        //! @brief name of the class
        virtual
        std::string
        name() const override { return "Materialized"+name_; }

        //! @brief number of stored function values
        auto
        datasize() const {return data_.size();}

    private:
        FiniteDomainVector domains_;
        data_t data_;
        fun_value_t zero_;
        std::string name_;

        template<class T>
        auto unique( const std::vector<T> &xs ) {
            auto vec = xs;
            std::sort( vec.begin(), vec.end() );
            vec.erase( std::unique( vec.begin(), vec.end() ), vec.end() );
            return vec;
        }

        auto
        extract_domains(const FiniteDomainVector &v) {
            auto ds = FiniteDomainVector();
            for ( auto x: this->vars() ) {
                ds.push_back(v[x]);
            }
            return ds;
        }


        //! @brief unique index calculated from the values of the
        //! function vars in the assignment
        auto
        index( const assignment_t & a ) const {
            size_t x = 0;
            size_t i = 0;
            auto vars = this->vars();
            for ( size_t i=0; i < domains_.size(); i++) {
                x *= domains_[i].size();
                x += a[vars[i]] - domains_[i].lb();
            }
            return x;
        }

        //! @brief size required for (non-sparse) vector of function
        //! values
        auto
        calc_size() const {
            return std::accumulate( domains_.begin(), domains_.end(), 1,
                       [](int acc, const FiniteDomain &x) { return x.size()*acc; } );
        }
    };

    //! A constraint is simply a boolean-valued function
    using Constraint = Function<bool>;
}

#endif
