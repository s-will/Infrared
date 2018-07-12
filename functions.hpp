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

#include <vector>
#include <memory>

#include "assignment.hpp"

namespace ired {

    /**
     * Dependencies specify a dependency between variables
     */
    class Dependency {
    public:
        using var_idx_t = int;

        Dependency(const std::vector<var_idx_t> &vars) : vars_(vars) {}

        const std::vector<var_idx_t> &
        vars() const {return vars_;}

    private:
        const std::vector<var_idx_t> vars_;
    };

    /**
     * Functions evaluate assignments of a subset of variables
     */
    template<class FunValue>
    class Function : public Dependency {
    public:
        using self_t = Function<FunValue>;
        using base_t = self_t;
        using assignment_t = Assignment;
        using fun_value_t = FunValue;

        Function(const std::vector<var_idx_t> &vars) : Dependency(vars) {}

        Function(const Function &fun) : Dependency(fun) {}

        virtual
        fun_value_t
        operator () (const assignment_t &) const {
            return fun_value_t();
        }

        virtual
        std::unique_ptr<base_t>
        clone() const {
            return std::make_unique<self_t>(*this);
        }

        virtual
        ~Function() {}
    };

    template<class FunValue, class UnderlyingFunction>
    class DynamicFunction : public Function<FunValue> {
    public:
        using self_t = DynamicFunction<FunValue,UnderlyingFunction>;
        using parent_t = Function<FunValue>;
        using base_t = typename parent_t::base_t;
        using var_idx_t = typename parent_t::var_idx_t;
        using assignment_t = typename parent_t::assignment_t;
        using fun_value_t = FunValue;


        template<class F>
        DynamicFunction(const std::vector<var_idx_t> &vars, const F &f) : parent_t(vars), f_(f) {}

        fun_value_t
        operator () (const assignment_t & a) const override { return f_(a); }

        std::unique_ptr<base_t>
        clone() const override {
            return std::make_unique<self_t>(*this);
        }

    private:
        using function_t = UnderlyingFunction;
        const function_t f_;
    };

    //! MaterializedFunction is used for the computed messages (DP
    //! matrices)
    //!
    //! @note via template parameter Container, one can
    //! choose to use a sparse data structure like
    //! Container=std::unordered_map<FunValue>
    //
    template<class FunValue, class Container=std::vector<FunValue> >
    class MaterializedFunction : public Function<FunValue> {
    public:
        using self_t = MaterializedFunction<FunValue,Container>;
        using parent_t = Function<FunValue>;
        using base_t = typename parent_t::base_t;

        using var_idx_t = typename parent_t::var_idx_t;
        using assignment_t = typename parent_t::assignment_t;
        using fun_value_t = FunValue;
        using data_t = Container;

        MaterializedFunction(const std::vector<var_idx_t> &vars, const std::vector<int> &domsizes) :
            parent_t(vars),
            domsizes_(domsizes),
            data_( calc_size() )
        {}

        fun_value_t
        operator () ( const assignment_t & a ) const override { return data_[ index_(a) ]; }

        std::unique_ptr<base_t>
        clone() const override {
            return std::make_unique<self_t>(*this);
        }

        void
        set( const assignment_t & a, const fun_value_t &val) {
            data_[ index_(a) ] = val;
        }

    private:

        data_t data_;

        const std::vector<int> domsizes_;

        int
        index_( const assignment_t & a ) const {
            int x = 0;

            for ( auto var : this->vars()) {
                x *= domsizes_[var];
                x += a[var];
            }

            return x;
        }

        auto
        calc_size() const {
            int x=1;

            for ( auto var : this->vars()) {
                x *= domsizes_[var];
            }
            return x;
        }
    };

    using Constraint = Function<bool>;

    template<class FunValue, class UnderlyingFunction>
    using DynamicConstraint = DynamicFunction<bool, UnderlyingFunction>;
}

#endif
