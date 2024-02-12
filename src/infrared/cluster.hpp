#ifndef INFRARED_CLUSTER_HPP
#define INFRARED_CLUSTER_HPP

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
 * @brief Defines the clusters (=bags) of the cluster tree.
 */


#include <vector>
#include <algorithm>

#include "functions.hpp"

namespace ired {
    /**
     * @brief Cluster (or bag) in the cluster tree
     *
     * A cluster consists of collections
     *  - variables
     *  - constraints
     *  - functions
     *
     * Implicitly, a cluster belongs to a cluster tree and, thus, a
     * constraint network CN. Variables are stored/identified by their unique
     * indices in the CN.
     */
    template<class FunValue=double>
    class Cluster {
    public:
        using var_idx_t = int;
        using function_t = Function<FunValue>;
        using constraint_t = Constraint;
        using assignment_t = Assignment;

        /**
         * @brief empty constructor
        */
        Cluster() {
        }

        /**
         * @brief constructor with variables
         *
         * @param vars vector of indices of variables in the cluster
         */
        explicit
        Cluster(const std::vector<var_idx_t> &vars)
            : vars_(vars) {
        }

        ~Cluster() {}

        /**
         * @brief Get variables
         *
         * @return vector of indices of the variable in the cluster
         */
        const std::vector<var_idx_t> &
        vars() const {return vars_;}

        /**
         * @brief Check whether cluster is empty
         */
        bool
        empty() const {
            return vars_.empty();
        }

        /**
         * @brief Get constraints
         *
         * @return vector of constraints assigned to the cluster
         */
        const auto &
        constraints() const {
            return constrs_;
        }

        /**
         * @brief Get functions
         *
         * @return vector of functions assigned to the cluster
         */
        const auto &
        functions() const {return funs_;}

        /**
         * @brief Assign constraint to this cluster
         */
        void
        add_constraint(const constraint_t *c) {
            constrs_.push_back(c);
        }

        /**
         * @brief Assign function to this cluster
         */
        void
        add_function(const function_t *f) {
            funs_.push_back(f);
        }

        void
        modify_function(size_t i, function_t *f) const {
            const_cast<std::vector<const function_t*>*>(&funs_)->operator[](i) = f;
        }

        /**
         * @brief Calculate separator variables
         *
         * @param parent parent cluster of this object
         * @return vector of the indices of the separator variables
         *
         * The separator variables are the variables of this object
         * that are already present in the parent. In other terms, the
         * cut set of the variables of this and the parent cluster
         */
        auto
        sep_vars(const Cluster &parent) const {
            auto sepvars = std::vector<var_idx_t>();
            const auto &pvars = parent.vars();
            for ( auto &var : vars_ ) {
                if (std::find(pvars.begin(), pvars.end(),var)
                    != pvars.end()) {
                    sepvars.push_back(var);
                }
            }
            return sepvars;
        }

        /**
         * @brief Calculate difference variables
         *
         * @param parent parent cluster of this object
         * @return vector of the indices of the difference variables
         *
         * The difference variables are the variables of this object that
         * are *not* already present in the parent.
         */
        auto
        diff_vars(const Cluster &parent) const {
            auto diffvars = std::vector<var_idx_t>();
            const auto &pvars = parent.vars();
            for ( auto &var : vars_ ) {
                if (std::find(pvars.begin(), pvars.end(),var)
                    == pvars.end()) {
                    diffvars.push_back(var);
                }
            }
            return diffvars;
        }

    private:
        std::vector<var_idx_t> vars_;
        std::vector<const constraint_t *> constrs_;
        std::vector<const function_t *> funs_;
    }; // end class Cluster
}

#endif
