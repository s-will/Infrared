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
#include <unordered_map>
#include <memory>

#include "assignment.hpp"

namespace ired {

    /**
     * Dependencies specify a dependency between variables
     *
     * @todo what does the engine need? Could we, e.g., rather offer
     * only a depends_on(var_idx_t) than accessor vars()?  Currently,
     * this almost requires functions to hold the vars as vector,
     * which e.g. makes function definitions more ugly.
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
    template<class FunValue=double>
    class Function : public Dependency {
    public:
        using self_t = Function<FunValue>;
        using base_t = self_t;
        using assignment_t = Assignment;
        using fun_value_t = FunValue;

        Function(const std::vector<var_idx_t> &vars) : Dependency(vars) {}

        // Function(const Function &fun) : Dependency(fun) {}

        virtual
        fun_value_t
        operator () (const assignment_t &) const = 0;

        virtual
        ~Function() {}
    };

    //! MaterializedFunction is used for the computed messages (DP
    //! matrices)
    //!
    //! @note via template parameter Container, one can
    //! choose to use a sparse data structure like
    //! Container=std::unordered_map<FunValue>
    //
    struct mapS {};
    struct vecS {};
    template<class FunValue, class T>
    struct container_selector {
        using type = T;
    };
    template<class FunValue>
    struct container_selector<FunValue,vecS> {
        using type = std::vector<FunValue>;
        static void init(type &x, size_t size, const FunValue &zero) {
            x.resize(size);
            std::fill(x.begin(),x.end(),zero);
        }
        static const FunValue& get(const type &x, size_t i, const FunValue &zero) {return x[i];}
        //static void set(type &x,size_t i, const FunValue &v) {x[i]=v;}
    };
    template<class FunValue>
    struct container_selector<FunValue,mapS> {
        using type = std::unordered_map<int,FunValue>;
        static void init(type &x, size_t size, const FunValue &zero) {
        }
        static const FunValue& get(const type &x, size_t i, const FunValue &zero) {
            auto it = x.find(i);
            if (it != x.end()) return it->second; else return zero;
        }
        //static void set(type &x,size_t i, const FunValue &v) {x[i]=v;}
    };

    template< class FunValue, class ContainerS=mapS >
    class MaterializedFunction : public Function<FunValue> {
    public:
        using self_t = MaterializedFunction<FunValue,ContainerS>;
        using parent_t = Function<FunValue>;
        using base_t = typename parent_t::base_t;

        using var_idx_t = typename parent_t::var_idx_t;
        using assignment_t = typename parent_t::assignment_t;
        using fun_value_t = FunValue;
        
        using data_t = typename container_selector<FunValue,ContainerS>::type;        

        MaterializedFunction(const std::vector<var_idx_t> &vars,
                             const std::vector<int> &domsizes, 
                             const fun_value_t &zero = fun_value_t()) :
            parent_t(vars),
            domsizes_(domsizes),
            data_(),
            zero_(zero)
        {
            container_selector<FunValue,ContainerS>::init(data_, calc_size(), zero_);
        }

        fun_value_t
        operator () ( const assignment_t & a ) const override {
            return container_selector<FunValue,ContainerS>::get(data_, index(a), zero_);
        }

        void
        set( const assignment_t & a, const fun_value_t &val) {
            data_[ index(a) ] = val;
        }

        auto
        datasize() const {return data_.size();}

    private:
        const std::vector<int> domsizes_;
        data_t data_;
        fun_value_t zero_;

        auto
        index( const assignment_t & a ) const {
            size_t x = 0;

            for ( auto var : this->vars()) {
                x *= domsizes_[var];
                x += a[var];
            }

            return x;
        }

        auto
        calc_size() const {
            size_t x=1;

            for ( auto var : this->vars()) {
                x *= domsizes_[var];
            }

            return x;
        }
    };

    using Constraint = Function<bool>;
}

#endif
