#ifndef INFRARED_ASSIGNMENT_HPP
#define INFRARED_ASSIGNMENT_HPP

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
#include <limits>
#include <cassert>
#include <iostream>

namespace ired {
    template<class EvaluationPolicy>
    class AssignmentIterator;

    // (partial) assignment assigns cn variables to values
    class Assignment {
    public:

        using assignment_t = Assignment;

        using var_idx_t = int;
        using var_value_t = int;
        static const int Undetermined = -1;

        //! @brief construct 'empty' with size
        Assignment(int size)
            :values_(size, Undetermined)
        {
        }

        decltype(auto)
        operator [](var_idx_t i) const {return values_[i];}

        decltype(auto)
        operator [](var_idx_t i) {return values_[i];}

        const auto &
        values() const { return values_; }

        //! @brief is a variable determined by the assignment?
        auto
        is_det(var_idx_t i) const {
            return values_[i] != Undetermined;
        }

        //! @brief set variables to undetermined
        void
        set_undet(const std::vector<var_idx_t> &xs) {
            for (auto x: xs) {
                values_[x] = Undetermined;
            }
        }

        /**
         * @brief make sub-assignment iterator
         * @see class AssignmentIterator
         *
         * During the iteration, constraints and functions are
         * evaluated.  We evaluate only those constraints/functions
         * that either have only dependencies on vars or have at least
         * on undetermined dependency on a variable in vars, where all
         * dependencies outside of vars are already determined.
         */
        template<class EP> // EP=EvaluationPolicy
        auto
        make_iterator(const std::vector<var_idx_t> &vars,
                      const std::vector<int> &domsizes,
                      const std::vector<const typename EP::constraint_t *> &
                      constraints,
                      const std::vector<const typename EP::function_t *> &
                      functions,
                      const EP &,
                      const typename EP::fun_value_t & initial_value = EP::one()
                      );


        /**
         * @brief evaluate constraints/functions where all variables
         * are already determined
         * @returns product of determined function values
         */
        template<class Function, class EP>
        auto
        eval_determined(const std::vector<const Function *> &functions, const EP &) {
            // for each function
            auto x = EP::one();
            for( const auto &f : functions ) {
                if ( all_of( f->vars().begin(), f->vars().end(), [&] (auto x) {return is_det(x);} ) )
                    {
                        // evaluate f
                        x = EP::multiplies( x, (*f)(*this) );
                    }
            }
            return x;
        }

    private:
        std::vector<var_value_t> values_;
    };  // end class Assignment

    /**
     * @brief Iterate over the assignments of a subset of variables
     *
     * @note the assignment iterator works destructively on parent assignment!
     *
     * The job of the assignment iterator goes beyond pure enumeration
     * of the sub-assignments that are valid due to the domains of a
     * subset of variables. It also checks constraints and evaluates
     * functions (the latter is still @todo).
     *
     * Constraints and functions are evaluated as soon as all their
     * variables are determined. For this purpose, the iterator
     * precomputes tables ('boards').
     * Function values are combined via EvaluationPolicy::multiplies; partial products
     * are tracked on the enumeration stack.
     */
    template<class EvaluationPolicy>
    class AssignmentIterator {
    public:
        using assignment_t = Assignment;
        using var_idx_t = assignment_t::var_idx_t;
        using var_value_t = assignment_t::var_value_t;
        using ep = EvaluationPolicy;
        using fun_value_t = typename ep::fun_value_t;
        using function_t = typename ep::function_t;
        using constraint_t = typename ep::constraint_t;

        AssignmentIterator(assignment_t &a,
                           const std::vector<var_idx_t> &vars,
                           const std::vector<int> &domsizes,
                           const std::vector<const constraint_t *> &constraints,
                           const std::vector<const function_t *> &functions,
                           const fun_value_t &initial_value
                           )
            : a_(a),
              vars_(vars),
              domsizes_(domsizes),
              top_( 0 )
        {
            for (auto var : vars_) {
                a_[var] = assignment_t::Undetermined;
            }

            // set up value stack
            value_stack_.resize(vars_.size()+1);
            value_stack_[0] = initial_value;

            if ( value_stack_[0]==ep::zero() ) {
                top_=-1;
                return;
            }


            constraint_board_ = create_board<constraint_t>(constraints);
            function_board_ = create_board<function_t>(functions);

            if ( vars_.size() > 0 ) {
                this -> operator ++();
            }
        }

        fun_value_t
        value() {
            assert(top_>=0);
            assert( ( top_==0 && vars_.size()==0 )
                    || top_ == static_cast<int>(vars_.size())-1 );
            return value_stack_[vars_.size()];
        }

        //! check constraints and value for current top_ entry
        bool
        is_valid_at_top() const {
            assert(top_ < (int)vars_.size());

            if ( value_stack_[top_+1] == ep::zero() )
                return false;

            //check all constraints at top_
            auto constraints = constraint_board_[top_];
            return all_of( constraints.begin(), constraints.end(),
                           [&](auto c) { return (*c)(a_); } );
        }

        //! compute function value due to assignment of top variable;
        //! put result on stack
        void
        eval_at_top() {
            assert(top_>=0);
            auto x = value_stack_[top_];
            for ( const auto f: function_board_[top_] ) {
                x = ep::multiplies( x, (*f)(a_) );
            }
            value_stack_[top_+1] = x;
        }

        /** @brief Next valid valuation
         *
         * Sets the assignment to the next valid valuation of
         * vars_
         *
         * Invariants:
         *   top_ points to the next to be enumerated variable in vars_
         *   top_==-1 means termination
         *   otherwise a_[vars_[top_]] is valid,
         *     if a_[vars_[top_]] < domsize_[vars_[top_]],
         *
         *   value_stack_[i] is the product of all function
         *   values, where the functions depend only on variables up
         *   to and including vars_[i-1], for 1<=i<=top_. If
         *   vars[top_] is undetermined, then value_stack_[top_+1] is undefined.
         *
         *   value_stack_[0] is always ep::one()
         */
        auto
        operator ++ () {
            //terminate if stack is empty
            if ( top_ < 0 ) {
                return *this;
            }
            while (top_ < static_cast<int>(vars_.size())) {
                assert(top_>=0);

                // determine the next valid candidate partial assignment
                do {
                    a_[vars_[top_]]++;

                    while ( a_[vars_[top_]] >= domsizes_[vars_[top_]] ) {
                        a_[vars_[top_]] = assignment_t::Undetermined;
                        top_ --;
                        if ( top_ < 0 ) {
                            // terminate if stack is empty
                            return *this;
                        }
                        a_[vars_[top_]] ++;
                    }

                    eval_at_top();

                    //repeat until the top variable has a valid value
                } while ( ! is_valid_at_top() );
                top_++;
            }
            top_--;
            assert( top_ == static_cast<int>(vars_.size())-1 );

            return *this;
        }

        /** @brief Check for termination
         */
        bool finished() {
            return top_ < 0;
        }

    private:
        template<class Function>
        struct board_t { using type = std::vector<std::vector<const Function *>>; };

        assignment_t &a_;
        const std::vector<var_idx_t> &vars_;
        const std::vector<int> &domsizes_;
        typename board_t<constraint_t>::type constraint_board_;
        typename board_t<function_t>::type function_board_;
        int top_;

        std::vector<fun_value_t> value_stack_;

        /**
         * @brief create the functions board
         *
         * The board is generated such that directly after fixing
         * variable i, one evaluates the functions in board[i].  It is
         * guaranteed that at this point, all arguments of the
         * functions are determined.  The board comprises only
         * functions, where at least one value in vars_ of the
         * iterator and vars of the function was undetermined before
         * start of the enumeration. (Recall that constraints are
         * boolean-valued functions).
         */
        template<class Function>
        auto
        create_board(const std::vector<const Function *> &functions) {

            typename board_t<Function>::type board(vars_.size());

            // for each function
            for( const auto &f : functions ) {

                // find 'last' undetermined variable of f
                // by iterating over the variables of f
                int last_idx = -1;
                for ( auto var : f->vars() ) {
                    // skip determined variables of f
                    if ( a_[var] != assignment_t::Undetermined ) continue;

                    // determine maximum index of this variable (in the order of the variables
                    // of the assignment iterator)
                    //    determine the index of the constraint variable in the iterator variables
                    auto idx = std::distance(vars_.begin(), find(vars_.begin(),vars_.end(),var));
                    // if the undetermined variable does not occur, then it will not be fixed during this enumeration --> we will not register f

                    //    determine maximum as last_idx
                    last_idx = std::max(last_idx, static_cast<int>(idx));
                }

                // here, ignore f if either it does not have undetermined variables, or
                // it has undetermined vars, but at least one is outside of the iterator vars
                if ( 0 <= last_idx && last_idx < static_cast<int>(vars_.size()) ) {
                    // register f at last_idx
                    board[last_idx].push_back( f );
                }
            }
            return board;
        }
    }; // end class AssignmentIterator

    const int Assignment::Undetermined;

    template<class EP>
    auto
    Assignment::make_iterator(const std::vector<var_idx_t> &vars,
                              const std::vector<int> &domsizes,
                              const std::vector<const typename EP::constraint_t *> &
                              constraints,
                              const std::vector<const typename EP::function_t *> &
                              functions,
                              const EP &,
                              const typename EP::fun_value_t & initial_value
                              ) {
        return AssignmentIterator<EP>(*this, vars, domsizes, constraints, functions, initial_value);
    }

}

#endif
