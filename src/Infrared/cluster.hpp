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


#include <vector>
#include <algorithm>

#include "functions.hpp"

namespace ired {
    /**
     * Cluster in the cluster tree decomposition
     *
     * A cluster consists of collections
     *  - variables
     *  - constraints
     *  - functions
     */
    template<class FunValue=double>
    class Cluster {
    public:
        using var_idx_t = int;
        using function_t = Function<FunValue>;
        using constraint_t = Constraint;
        using assignment_t = Assignment;


        Cluster() {
        }

        Cluster(const std::vector<var_idx_t> &vars)
            : vars_(vars) {
        }

        Cluster(const std::vector<var_idx_t> &vars,
                const std::vector<constraint_t *> &constraints,
                const std::vector<function_t *> &functions
                )
            : vars_(vars),
              constrs_(constraints),
              funs_(functions) {
        }

        const std::vector<var_idx_t> &
        vars() const {return vars_;}

        const auto &
        constraints() const {
            return constrs_;
        }

        const auto &
        functions() const {return funs_;}

        void
        add_variable(const var_idx_t &x) {
            vars_.push_back(x);
        }

        void
        add_constraint(const constraint_t *c) {
            constrs_.push_back(c);
        }

        void
        add_function(const function_t *f) {
            funs_.push_back(f);
        }

        auto
        sep_vars(const Cluster &parent) const {
            auto vars = std::vector<var_idx_t>();
            const auto &pvars = parent.vars();
            for ( auto &var : vars_ ) {
                if (std::find(pvars.begin(), pvars.end(),var)
                    != pvars.end()) {
                    vars.push_back(var);
                }
            }
            return vars;
        }

        auto
        diff_vars(const Cluster &parent) const {
            auto vars = std::vector<var_idx_t>();
            const auto &pvars = parent.vars();
            for ( auto &var : vars_ ) {
                if (std::find(pvars.begin(), pvars.end(),var)
                    == pvars.end()) {
                    vars.push_back(var);
                }
            }
            return vars;
        }

    private:
        std::vector<var_idx_t> vars_;
        std::vector<const constraint_t *> constrs_;
        std::vector<const function_t *> funs_;
    };
}

#endif
