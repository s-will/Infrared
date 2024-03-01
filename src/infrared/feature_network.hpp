#ifndef INFRARED_FEATURE_NETWORK_HPP
#define INFRARED_FEATURE_NETWORK_HPP

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
#include <limits>
#include <iterator>
#include <cstdlib>

#include "cluster.hpp"
#include "finite_domain.hpp"

namespace ired {

    /**
     * @brief The evaluation policy to calculate partition
     * functions.
     *
     * An evaluation policy defines an algebra with operations
     * mul and plus, as well as respective neutral elements one
     * and zero. The following must hold
     *
     * mul is associative and symmetric
     * plus is associative and symmetric
     * mul(one,x) == x
     * mul(zero,x) == 0.0
     * plus(zero,x) == x
     *
     * In the standard policy mul corresponds to *, plus to +,
     * zero to 0.0, and one to 1.0
     */
    template<class FunValue>
    class PFEvaluationPolicy {
    public:
        using fun_value_t = FunValue;
        using constraint_t = Function<bool>;
        using function_t = Function<fun_value_t>;

        // selector for stochastic backtracking
        class selector {
        public:
            selector(const fun_value_t &value): 
                value_(value),
                r_( rand()/(RAND_MAX+1.0) * value) {
            }

            auto
            select(const fun_value_t &x) {
                return x > r_;
            }
        private:
            fun_value_t value_;
            fun_value_t r_;
        };

        static
        fun_value_t
        plus(const fun_value_t &x, const fun_value_t &y) {
            return x+y;
        }

        static
        fun_value_t
        minus(const fun_value_t &x, const fun_value_t &y) {
            return x-y;
        }

        static
        fun_value_t
        mul(const fun_value_t &x, const fun_value_t &y) {
            return x*y;
        }

        static
        fun_value_t
        del(const fun_value_t &x, const fun_value_t &y) {
            return x/y;
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

        static
        fun_value_t
        almost_zero(const fun_value_t &x) {
            return std::abs(x)<1e-14; // TO DO: adjust the exact value
        }
    };

    //! define PFEvaluationPolicy as standard
    template<class FunValue>
    class StdEvaluationPolicy: public PFEvaluationPolicy<FunValue> {
    };

    /**
     * @brief Evaluation Strategy for Optimization (max/+); defining the arctic semiring
     *
     * @see StdEvaluationPolicy
     *
     * Defines the arctic semiring
     *
     * mul corresponds to +, plus to max,
     * zero to -infty, and one to 0
     *
     */
    template<class FunValue>
    class ArcticEvaluationPolicy {
    public:
        using fun_value_t = FunValue;
        using constraint_t = Function<bool>;
        using function_t = Function<fun_value_t>;

        // selector for 'classic' backtracking
        class selector {
        public:
            selector(const fun_value_t &value): value_(value) {}
            auto
            select(const fun_value_t &x) {
                return x==value_;
            }
        private:
            fun_value_t value_;
        };

        static
        fun_value_t
        plus(const fun_value_t &x, const fun_value_t &y) {
            return std::max(x,y);
        }

        static
        fun_value_t
        mul(const fun_value_t &x, const fun_value_t &y) {
            return x + y;
        }

        static
        fun_value_t
        one() {
            return fun_value_t();
        }

        static
        fun_value_t
        zero() {
            return std::numeric_limits<fun_value_t>::min();
        }
    };

    /**
     * @brief The evaluation policy to combine constraint values.  Here,
     * mul corresponds to &&, plus to ||, one to true, and zero
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
        mul(const fun_value_t &x, const fun_value_t &y) {
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
     * @brief the feature network
     *
     * A feature network consists of sets of
     *  - variables
     *  - constraints
     *  - functions
     *
     * The variables are indexed (0..n-1). Each variable has a finite
     * domain.  An object holds shared pointers to the functions and
     * constraints of the network. In this way, it guarantees their
     * existence as long as the network exists.
     */
    template<class FunValue, class EvaluationPolicy> class FeatureNetwork {
    public: using var_idx_t = int; using var_value_t = int;

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
        FeatureNetwork() {};

        /**
         * @brief Construct with domains
         */
        explicit
        FeatureNetwork(const FiniteDomainVector &domains)
            : domains_(domains) {
        };

        /**
         * @brief Construct with domains of equal size
         */
        FeatureNetwork(int num_vars, FiniteDomain domain)
            : domains_(num_vars, domain) {
        };


        ~FeatureNetwork() {};


        /**
         * @brief add variable with domain size
         *
         * @returns 0-based index of the new variable
         */
        var_idx_t
        add_variable(FiniteDomain domain) {
            domains_.push_back(domain);
            return domains_.size()-1;
        }

        /**
         * @brief add constraint
         *
         * @param x shared pointer to constraint
         * @return the added constraint
         *
         * If the virtual method x->auto_materialize() returns true,
         * then the constraint *x is materialized as mx and the cn holds
         * a shared pointer to *mx. Otherwise, the network holds a
         * shared pointer to *x.
         */
        auto
        add_constraint(const std::shared_ptr<constraint_t> &x) {
            if ( x->auto_materialize() ) {
                //materialize it

                // for evaluating the constraints, construct a
                // constraint network with the same variables
                // and domains as this one, which evaluates constraints
                auto cn = FeatureNetwork<bool, StdEvaluationPolicy<bool>>(domains_);
                auto mx = std::make_shared<MaterializedFunction<bool, vecS>>(x.get(),cn);

                constraints_.push_back( mx );
            } else {
                constraints_.push_back( x );
            }
            return constraints_.back().get();
        }

        /**
         * @brief add function
         *
         * @param x shared pointer to function
         * @return the added function
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
        num_vars() const {return domains_.size();}

        /**
         * @brief Get vector of domains (read only)
         */
        const
        auto & domains() const {
            return domains_;
        }

    private:
        FiniteDomainVector domains_;

        std::vector<std::shared_ptr<function_t>> functions_;
        std::vector<std::shared_ptr<constraint_t>> constraints_;
    };

}
#endif
