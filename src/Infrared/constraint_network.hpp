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


#include <vector>
#include <memory>

#include "cluster.hpp"

namespace ired {
    /** A constraint network consists of sets of
     *  - variables
     *  - constraints
     *  - functions
     *
     * @note The variables in a contraint network are indexed (0..n-1). Each variable has a domain of values 0..domsize
     */
    template<class FunValue=double>
    class ConstraintNetwork {
    public:
        using var_idx_t = int;
        using var_value_t = int;

        using fun_value_t = FunValue;

        using assignment_t = Assignment;
        using function_t = Function<FunValue>;
        using constraint_t = Constraint;

        //! @brief suitable cluster type for the constraint network
        using cluster_t = Cluster<fun_value_t>;

        /**
         * @brief Construct empty
         */
        ConstraintNetwork() {};

        /**
         * @brief Construct with domains
         */
        ConstraintNetwork(std::vector<int> &domsizes)
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

        // /**
        //  * @brief add constraint (cloning)
        //  * @returns pointer to the new constraint
        //  */
        // auto
        // add_constraint(const constraint_t &x) {
        //     return add_constraint( x.clone() );
        // }

        // /**
        //  * @brief add function
        //  */
        // auto
        // add_function(const function_t &x)  {
        //     return add_function( x.clone() );
        // }

        /**
         * @brief add constraint (moving ownership)
         */
        auto
        add_constraint(const std::shared_ptr<constraint_t> &x) {
            constraints_.push_back( x );
            return constraints_.back().get();
        }

        /**
         * @brief add function
         */
        auto
        add_function(const std::shared_ptr<function_t> &x) {
            functions_.push_back( x );
            return functions_.back().get();
        }


        int
        num_vars() const {return domsizes_.size();}

        const auto & domsizes() const {
            return domsizes_;
        }

    private:
        std::vector<int> domsizes_;

        std::vector<std::shared_ptr<function_t>> functions_;
        std::vector<std::shared_ptr<constraint_t>> constraints_;


    };


    //NOTE: it should be possible to directly (manually) construct a tree decomposition
    // such that the cn is filled implicitely

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

}
#endif
