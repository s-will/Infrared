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

/**
 * @file
 *
 * @brief Defines assignments of the variables of a constraint network to values
 */

#include <vector>
#include <limits>
#include <cassert>
#include <iostream>
#include <functional>


#include "finite_domain.hpp"

namespace ired {
    template<class FeatureNetwork>
    class AssignmentIterator;

    //@brief Vector type selector "no bit vector"
    //!
    //! define a selector for std::vector while avoiding the bit-vector
    //! implementation for vector<bool>
    //!
    //! This was introduced to fix a bug: vector<bool> as value stack does not work for yet
    //! unknown reason; the use of vector<int> resolves the problem
    template<class T>
    struct vector_nbv_sel {
        using type = std::vector<T>;
    };

    template<>
    struct vector_nbv_sel<bool> {
        using type = std::vector<char>;
    };


    /** @brief A (partial) assignment of variables to values
     *
     * An object of Assignment associates values to a series of
     * /finite domain/ variables (indexed 0..n-1); typically the
     * variables belong to a constraint network CN (@see
     * FeatureNetwork).
     *
     * In an assignment for the CN, variables are associated to values
     * in their finite integer domain (according to CN) or are
     * undetermined (encoded domain specifically!).
     * Assignments are total, if no variable is undetermined.
     *
     * Partial assignments can be enumerated using objects of class
     * AssignmentIterator. Such iterators are generated for an
     * assignment using make_iterator; typically they are used to
     * iterate over a subset of the variables in the CN while
     * evaluating functions and checking constraints.
     *
     * An assignment is valid for a set of constraints, iff all the
     * constraints are satisfied (evaluate to true).  The value of an
     * assignment for a set of functions is the product (as defined by
     * the EvaluationPolicy) of the function values for the
     * assignment. Only assignments with non-zero (zero def'd by
     * EvaluationPolicy) evaluation are valid for the functions.
     */
    class Assignment {
    public:

        //! assignment type
        using assignment_t = Assignment;

        //! variable index type
        using var_idx_t = int;

        //! variable value type
        using var_value_t = int;

    public:
        /**
         * @brief construct empty
         *
         * @param domains domains of the variables; attention: domains must
         * stay available
         *
         * after construction, all variables are undetermined
         */
        explicit
        Assignment(const FiniteDomainVector &domains)
            :values_(domains.size()), domains_(domains)
        {
            for ( int i=0; i < values_.size(); ++i ) {
                set_undet(i);
            }
        }

        /**
         * @brief get value of a variable (read only)
         * @param i variable index
         * @return value of variable i
         */
        const auto &
        operator [](var_idx_t i) const {return values_[i];}

        /**
         * @brief get value of a variable (read/write)
         * @param i variable index
         * @return value of variable i
         */
        auto &
        operator [](var_idx_t i) {return values_[i];}

        /**
         * @brief get values (read only)
         */
        const auto &
        values() const { return values_; }

        /**
         * @brief get number of variables in assignment
         */
        auto
        size() const { return values_.size(); }

        /**
         * @brief check whether variable is determined
         * @param i variable index
         */
        auto
        is_det(var_idx_t i) const {
            return values_[i] != domains_[i].undet();
        }
 
        /**
         * @brief set variables to undetermined
         * @param i index of variable
         */
        void
        set_undet(var_idx_t i) {
            values_[i] = domains_[i].undet();
        }
 
        /**
         * @brief set variables to undetermined
         * @param xs vector of variable indices
         */
        void
        set_undet(const std::vector<var_idx_t> &xs) {
            for (auto x: xs) {
                set_undet(x);
            }
        }

        /**
         * @brief increment value of variable
         * @param i index of variable
         */
        void
        inc(var_idx_t i) {
            domains_[i].inc( values_[i] );
        }

        /**
         * @brief make sub-assignment iterator
         * @see class AssignmentIterator
         *
         * @param vars vector of indices of the to-be-enumerated variables
         * @param cn constraint network, where everything takes place
         * @param constraints check assignments
         * @param functions check and evaluate assignments
         */
        template<class FeatureNetwork>
        auto
        make_iterator(const std::vector<var_idx_t> &vars,
                      const FeatureNetwork &cn,
                      const std::vector<const typename FeatureNetwork::constraint_t *> &
                      constraints,
                      const std::vector<const typename FeatureNetwork::function_t *> &
                      functions,
                      const typename FeatureNetwork::fun_value_t & initial_value
                      );

        /**
         * @brief evaluate constraints/functions where all variables
         * are already determined
         *
         * @returns product of determined function values
         */
        template<class Function, class EP>
        auto
        eval_determined(const std::vector<const Function *> &functions, const EP &) {
            // for each function
            auto x = EP::one();
            for( const auto &f : functions ) {
                if ( all_of( f->vars().begin(), f->vars().end(),
                             [&] (auto v) { return is_det(v); } ) )
                    {
                        // evaluate f
                        x = EP::mul( x, (*f)(*this) );
                    }
            }
            return x;
        }

    private:
        /**
         * @brief vector holding a value (or Undetermined) for each
         * variable
         */
        std::vector<var_value_t> values_;
        const FiniteDomainVector &domains_;
    }; // end class Assignment



    /**
     * @brief Iterate over the assignments of a subset of variables
     *
     * @note the assignment iterator works destructively on parent assignment!
     *
     * The job of the assignment iterator goes beyond pure enumeration
     * of the sub-assignments that are valid due to the domains of a
     * subset of variables. It also checks validity due to constraints
     * and functions and evaluates the functions.
     *
     * Constraints and functions are evaluated as soon as all their
     * variables are determined. For this purpose, the iterator
     * precomputes tables ('boards').  Function values are combined
     * via EvaluationPolicy::mul; partial products are tracked
     * on the enumeration stack.
     */
    template<class FeatureNetwork>
    class AssignmentIterator {
    public:
        using assignment_t = Assignment;
        using var_idx_t = assignment_t::var_idx_t;
        using var_value_t = assignment_t::var_value_t;
        using ep = typename FeatureNetwork::evaluation_policy_t;
        using fun_value_t = typename ep::fun_value_t;
        using function_t = typename ep::function_t;
        using constraint_t = typename ep::constraint_t;

        /**
         * @brief constructor
         *
         * @param vars vector of indices of the to-be-enumerated variables
         * @param cn constraint network, where everything takes place
         * @param constraints check assignments
         * @param functions check and evaluate assignments
         *
         * During the iteration, constraints and functions are
         * evaluated.  We evaluate only those constraints/functions
         * that have at least on undetermined dependency on a variable
         * in vars, where all dependencies outside of vars are already
         * determined.
         *
         * After construction with non-empty vars, the assignment is
         * either set to the first valid assignment or finished() is
         * true. After construction with empty vars, assignment is
         * unchanged and finished() is false, but will be true
         * immediately after a call to operator ++()!
         */
        AssignmentIterator(assignment_t &a,
                           const std::vector<var_idx_t> &vars,
                           const FeatureNetwork &cn,
                           const std::vector<const constraint_t *> &constraints,
                           const std::vector<const function_t *> &functions,
                           const fun_value_t &initial_value
                           )
            : a_(a),
              vars_(vars),
              domains_(cn.domains()),
              initial_value_(initial_value),
              top_(0),
              stage1_size_(-1) // init without staging
        {
            a_.set_undet( vars_ );

            constraint_board_ = create_board<constraint_t>(constraints);
            function_board_ = create_board<function_t>(functions);

            this->reset();
        }

        /**
         * @brief Reset iterator
         */
        void reset() {
            a_.set_undet( vars_ );

            // set up value stack
            value_stack_.resize(vars_.size()+1);

            value_stack_[0] = initial_value_;

            if ( value_stack_[0]==ep::zero() ) {
                top_=-1;
                return;
            }

            if ( vars_.size() > 0 ) {
                this -> operator ++();
            }
        }

        /**
         * @brief get evaluation according to functions
         *
         * after initialization or a call of operator++, this yields
         * the product of all functions on enumerated variables, where
         * at least one variable is determined only by the enumeration
         */
        auto
        value() {
            assert(top_>=0);
            assert( ( top_==0 && vars_.size()==0 )
                    || top_ == static_cast<int>(vars_.size())-1 );
            return value_stack_[vars_.size()];
        }

        /**
         * @brief register hook to run after finishing stage 2 enumeration
         * @param stage1_size number of variables in stage 1
         * @param hook function to be called after finishing stage 2 enumeration
         *
         * This is used in two-stage enumeration over two disjoint
         * sets of variables; in the evaluation of a CN, this allows
         * to marginalize over the second set and store the results as
         * value of a function (the 'message') of the variables in the
         * first set.
         */
        void
        register_finish_stage2_hook(int stage1_size, const std::function<void()> &hook) {
            stage1_size_ = stage1_size;
            finish_stage2_hook_ = hook;
        }


        /** @brief Next valid partial assignment
         *
         * Sets the assignment to the next valid assignment of vars_
         */
        const auto &
        operator ++() {
            /*
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

            //terminate if stack is empty
            if ( top_ < 0 ) {
                return *this;
            }
            while (top_ < static_cast<int>(vars_.size())) {
                assert(top_>=0);

                // determine the next valid candidate partial assignment
                do {
                    a_.inc( vars_[top_] );

                    while ( a_[vars_[top_]] > domains_[vars_[top_]].ub() ) {

                        // call hook after each complete enum of stage2
                        // (i.e. at return to stage1)
                        if (top_ == stage1_size_) {
                            finish_stage2_hook_();
                        }

                        a_.set_undet( vars_[top_] );

                        top_ --;
                        if ( top_ < 0 ) {
                            // terminate if stack is empty
                            return *this;
                        }
                        a_.inc( vars_[top_] );
                    }

                    //repeat until the top variable has a valid value

                     if ( ! constraints_valid_at_top() )
                         continue;

                    if ( ! functions_valid_at_top() )
                        continue;

                     eval_at_top();
                     if ( value_stack_[top_+1] == ep::zero() )
                         continue;

                     break;
                } while (true);
                top_++;
            }

            // this test is required for the case where diff.size()==0
            if (top_ == stage1_size_) {
                finish_stage2_hook_();
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
        const FiniteDomainVector &domains_;
        const fun_value_t initial_value_;
        typename board_t<constraint_t>::type constraint_board_;
        typename board_t<function_t>::type function_board_;
        int top_;


        using value_stack_t = typename vector_nbv_sel<fun_value_t>::type;
        value_stack_t value_stack_;

        int stage1_size_;
        std::function<void()> finish_stage2_hook_;


        //! @brief check constraints for current top_ entry
        bool
        constraints_valid_at_top() const {
            assert(top_ < (int)vars_.size());

            //check all constraints at top_
            const auto &constraints = constraint_board_[top_];
            return all_of( constraints.begin(), constraints.end(),
                           [&](auto c) { return (*c)(a_); } );
        }

        //! @brief check functions for current top_ entry
        //!
        //! @note here, 'valid' is used in the sloppy sense that no
        //! function is 'guaranteed' to have zero value (@see
        //! Function<>::guaranteed_zero())
        bool
        functions_valid_at_top() const {
            assert(top_ < (int)vars_.size());

            //check all constraints at top_
            const auto &functions = function_board_[top_];
            return all_of( functions.begin(), functions.end(),
                           [&](auto f) { return ! f->guaranteed_zero(a_); } );
        }

        //! compute function value due to assignment of top variable;
        //! put result on stack
        void
        eval_at_top() {
            assert(top_>=0);
            auto x = value_stack_[top_];
            for ( const auto f: function_board_[top_] ) {
                x = ep::mul( x, (*f)(a_) );
            }

            value_stack_[top_+1] = x;
        }

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
                    if ( a_.is_det(var) ) continue;

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

    // const int Assignment::Undetermined;

    template<class CN>
    auto
    Assignment::make_iterator(const std::vector<var_idx_t> &vars,
                              const CN &cn,
                              const std::vector<const typename CN::constraint_t *> &
                              constraints,
                              const std::vector<const typename CN::function_t *> &
                              functions,
                              const typename CN::fun_value_t & initial_value
                              ) {
        return AssignmentIterator<CN>(*this,
                                      vars,
                                      cn,
                                      constraints,
                                      functions,
                                      initial_value);
    }

}

#endif
