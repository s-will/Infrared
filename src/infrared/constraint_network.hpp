#ifndef INFRARED_CONSTRAINT_NETWORK_HPP
#define INFRARED_CONSTRAINT_NETWORK_HPP

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
 * @brief Defines constraint networks and standard evaluation
 * policies.
 */

#include <vector>
#include <memory>

#include "cluster.hpp"

namespace ired {

    /**
     * @brief The standard evaluation policy to calculate partition
     * functions.
     *
     * An evaluation policy defines an algebra with operations
     * multiplies and plus, as well as respective neutral elements one
     * and zero. The following must hold
     *
     * multiplies is associative and symmetric
     * plus is associative and symmetric
     * multiplies(one,x) == x
     * multiplies(zero,x) == x
     * plus(zero,x) == x
     *
     * In the standard policy multiplies corresponds to *, plus to +,
     * zero to 0.0, and one to 1.0
     */
    template<class FunValue>
    class StdEvaluationPolicy {
    public:
        using fun_value_t = FunValue;
        using constraint_t = Function<bool>;
        using function_t = Function<fun_value_t>;

        static
        fun_value_t
        plus(const fun_value_t &x, const fun_value_t &y) {
            return x+y;
        }

        static
        fun_value_t
        multiplies(const fun_value_t &x, const fun_value_t &y) {
            return x*y;
        }

        static
        fun_value_t
        one() {
            return fun_value_t(1);
        }

        static
        fun_value_t
        zero() {
            return fun_value_t();
        }
    };

    /**
     * @brief The evaluation policy to combine constraint values.  Here,
     * multiplies corresponds to &&, plus to ||, one to true, and zero
     * to false.
     */
    template<>
    class StdEvaluationPolicy<bool> {
    public:
        using fun_value_t = bool;
        using constraint_t = Constraint;
        using function_t = Constraint;

        static
        fun_value_t
        plus(const fun_value_t &x, const fun_value_t &y) {
            return x || y;
        }

        static
        fun_value_t
        multiplies(const fun_value_t &x, const fun_value_t &y) {
            return x && y;
        }

        static
        fun_value_t
        one() {
            return true;
        }

        static
        fun_value_t
        zero() {
            return false;
        }
    };

    /**
     * @brief the constraint network
     *
     * A constraint network consists of sets of
     *  - variables
     *  - constraints
     *  - functions
     *
     * The variables in a contraint network are indexed (0..n-1). Each
     * variable has a domain of values 0..domsize.  A constraint
     * network object holds shared pointers to the functions and
     * constraints of the network. In this way, it guarantees their
     * existence as long as the network exists.
     */
    template<class FunValue=double, class EvaluationPolicy=StdEvaluationPolicy<FunValue>>
    class ConstraintNetwork {
    public:
        using var_idx_t = int;
        using var_value_t = int;

        using fun_value_t = FunValue;

        using assignment_t = Assignment;
        using function_t = Function<FunValue>;
        using constraint_t = Constraint;

        using evaluation_policy_t = EvaluationPolicy;

        //! @brief suitable cluster type for the constraint network
        using cluster_t = Cluster<fun_value_t>;

        /**
         * @brief Construct empty
         */
        ConstraintNetwork() {};

        /**
         * @brief Construct with domains
         */
        ConstraintNetwork(const std::vector<int> &domsizes)
            : domsizes_(domsizes) {
        };

        /**
         * @brief Construct with uniform domains
         */
        ConstraintNetwork(int num_vars, int domsize)
            : domsizes_(num_vars, domsize) {
        };

        /**
         * @brief add variable with domain size
         *
         * @returns 0-based index of the new variable
         */
        var_idx_t
        add_variable(int domainsize) {
            domsizes_.push_back(domainsize);
            return domsizes_.size()-1;
        }

        /**
         * @brief add constraint
         *
         * @param x shared pointer to constraint
         *
         * The network holds a shared pointer to *x.
         */
        auto
        add_constraint(const std::shared_ptr<constraint_t> &x) {
            constraints_.push_back( x );
            return constraints_.back().get();
        }

        /**
         * @brief add function
         *
         * @param x shared pointer to function
         *
         * If the virtual method x->auto_materialize() returns true,
         * then the function *x is materialized as mx and the cn holds
         * a shared pointer to *mx. Otherwise, the network holds a
         * shared pointer to *x.
         */
        auto
        add_function(const std::shared_ptr<function_t> &x) {
            if ( x->auto_materialize() ) {
                //materialize it
                auto mx = std::make_shared<MaterializedFunction<typename function_t::fun_value_t, vecS>>(x.get(),*this);
                functions_.push_back( mx );
            } else {
                functions_.push_back( x );
            }
            return functions_.back().get();
        }

        /**
         * @brief Number of variables
         */
        int
        num_vars() const {return domsizes_.size();}

        /**
         * @brief Get vector of domain sizes (read only)
         */
        const
        auto & domsizes() const {
            return domsizes_;
        }

    private:
        std::vector<int> domsizes_;

        std::vector<std::shared_ptr<function_t>> functions_;
        std::vector<std::shared_ptr<constraint_t>> constraints_;


    };

}
#endif
