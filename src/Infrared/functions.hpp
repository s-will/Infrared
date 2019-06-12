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

#include "assignment.hpp"

namespace ired {

    template<class FunValue, class EP>
    class ConstraintNetwork;

    /**
     * @brief Dependencies specify a dependency between variables
     */
    class Dependency {
    public:
	//! type of variable index
        using var_idx_t = int;

	/**
	 * Constructor
	 * @param vars vector of variable indices
	 */
        Dependency(const std::vector<var_idx_t> &vars) : vars_(vars) {}

	/**
	 * @brief Return variable index vector
	 */
        const std::vector<var_idx_t> &
        vars() const {return vars_;}

    private:
        const std::vector<var_idx_t> vars_;
    };

    /**
     * @brief Functions evaluate assignments of a subset of variables
     * @tparam FunValue type of Function value. Default: double
     *
     * @note
     * The class with boolean FunValue is called Constraint
     * @code{.hpp}
     * using Constraint = Function<bool>
     * @endcode
     */
    template<class FunValue=double>
    class Function : public Dependency {
    public:
        using self_t = Function<FunValue>;
        using base_t = self_t;
        using assignment_t = Assignment;
        using fun_value_t = FunValue;

        Function(const std::vector<var_idx_t> &vars) : Dependency(vars) {}

        // Function(const Function &fun) : Dependency(fun) {}

	/** Function operator evaluating given assignment
	 *  @return evaluation
         */
        virtual
        fun_value_t
        operator () (const assignment_t &) const = 0;

        /** quick check whether function is definitely zero at assignment
         * @param a assignment
         * @return whether defined
         *
         *  some functions know how to quick check for zero (in particular sparse materialized functions)
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
     * @brief map selector class
     *
     * for selecting a sparse table to represent
     * materialized functions
     */
    struct mapS {};

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
        using type = std::vector<FunValue>;
        static void init(type &x, size_t size, const FunValue &zero) {
            x.resize(size);
            std::fill(x.begin(),x.end(),zero);
        }
        static bool guaranteed_zero(const type &x, size_t i) {
            return false;
        }
        static const FunValue& get(const type &x, size_t i, const FunValue &zero) {return x[i];}
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
     * @brief A materialized function
     *
     * A function that holds a table of function values. It is used to
     * represent computed messages (DP matrices) during evalution of
     * the cluster tree, or precomputed functions.
     *
     * Supports sparse and non-sparse representation of the function
     * value table. The choice is made by a template mechanism using
     * selector classes.
     */
    template< class FunValue, class ContainerS=mapS >
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
         * @param domsizes domain sizes of the variables in the cn; domsizes must outlife this object
         * @param zero value (default for not explicitly set function values)
         *
         * The (non-zero) function values are typically set after
         * construction, using method @see set()
         */
        template<class ConstraintNetwork>
        MaterializedFunction(const std::vector<var_idx_t> &vars,
                             const ConstraintNetwork &cn
                             )
            :
            parent_t(vars),
            domsizes_(cn.domsizes()),
            zero_(ConstraintNetwork::evaluation_policy_t::zero())
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
        template<class ConstraintNetwork>
        MaterializedFunction(const function_t *function,
                             const ConstraintNetwork &cn
                             )
            :
            parent_t(function->vars()),
            domsizes_(cn.domsizes()),
            zero_(ConstraintNetwork::evaluation_policy_t::zero())
        {
            using ep = typename ConstraintNetwork::evaluation_policy_t;

            container_selector<FunValue,ContainerS>::init(data_, calc_size(), zero_);

            auto a = Assignment(cn.num_vars());

            auto constraints = std::vector<const typename ConstraintNetwork::constraint_t *>();
            auto functions = std::vector<const typename ConstraintNetwork::function_t *> {function};

            auto it = a.make_iterator(function->vars(),
                                      cn,
                                      constraints,
                                      functions,
                                      a.eval_determined(functions, ep()) //evaluate 0-ary functions
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

        //! @brief whether to automatically materialzize when added to CN
        virtual
        bool
        auto_materialize() const override { return false; }

        //! @brief name of the class
        virtual
        std::string
        name() const override { return "MaterializedFunction"; }

        //! @brief number of stored function values
        auto
        datasize() const {return data_.size();}

    private:
        const std::vector<int> &domsizes_;
        data_t data_;
        fun_value_t zero_;

        //! @brief unique index calculated from the values of the
        //! function vars in the assignment
        auto
        index( const assignment_t & a ) const {
            size_t x = 0;
            for ( auto var : this->vars()) {
                x *= domsizes_[var];
                x += a[var];
            }
            return x;
        }

        //! @brief size required for (non-sparse) vector of function
        //! values
        auto
        calc_size() const {
            size_t x=1;
            for ( auto var : this->vars()) {
                x *= domsizes_[var];
            }
            return x;
        }
    };

    //! A constraint is simply a boolean-valued function
    using Constraint = Function<bool>;
}

#endif
