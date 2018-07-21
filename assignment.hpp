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
        is_det(var_idx_t i) {return values_[i] != Undetermined;}

        void
        reset(std::vector<var_idx_t> &vars) {
            for (auto var : vars) {
                values_[var] = Undetermined;
            }
        }

        // Iterate over possible assignments
        // assignment iterator works destructively on parent assignment!
        //
        // Special case: empty variable set: exactly one valid assignment
        template<class Constraint>
        class AssignmentIterator {
        public:
            using constraint_t = Constraint;

            AssignmentIterator(assignment_t &a,
                               const std::vector<var_idx_t> &vars,
                               const std::vector<int> &domsizes,
                               const std::vector<const Constraint *> &constraints)
                : a_(a),
                  vars_(vars),
                  domsizes_(domsizes),
                  top_( 0 )
            {
                for (auto var : vars_) {
                    a_[var] = Undetermined;
                }
                constraint_board_ = create_constraint_board(constraints);
                
                if ( vars_.size() > 0 ) {
                    this -> operator ++();
                }
            }

            // auto &
            // operator * () {
            //     return a_;
            // }


            // /**
            //  * @brief Reset iterator
            //  *
            //  * Resets the assignment such that the iteration starts again
            //  *
            //  * After reset, the valuation is not necessarily valid (use ++ operator!)
            //  */
            // auto
            // reset() {
            //     // the entries for vars_ in a are (ab)used as a stack
            //     // top_ points to the next to be determined variable
            //     // var=vars_[top_]; a_[var] is the next candidate
            //     // value
            //     for (auto var : vars_) {
            //         a_[var] = Undetermined;
            //     }
            //     top_ = std::min(static_cast<int>(vars_.size())-1,0);
            //     return *this;
            // }

            //! check constraints for current top_ entry
            bool
            is_valid() const {
                //check all constraints at top_
                auto constraints = constraint_board_[top_];

                return all_of( constraints.begin(), constraints.end(),
                               [&](auto c) { return (*c)(a_); } );
            }

            /** @brief Next valid valuation
             *
             * Sets the assignemnt to the next valid valuation of
             * vars_
             *
             * Invariants:
             *   top_ points to the next to be enumerated variable
             *   top_==-1 means termination
             *   otherwise a_[vars_[top_]] is valid, 
             *     if a_[vars_[top_]] < domsize_[vars_[top_]],
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
                            a_[vars_[top_]] = Undetermined;
                            top_ --;
                            if ( top_ < 0 ) {
                                // terminate if stack is empty
                                return *this;
                            }
                            a_[vars_[top_]] ++;
                        }
                        //repeat until is_valid
                    } while ( ! is_valid() );
                    
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

            using constraint_board_t = std::vector<std::vector<const constraint_t *>>;

            assignment_t &a_;
            const std::vector<var_idx_t> &vars_;
            const std::vector<int> &domsizes_;
            constraint_board_t constraint_board_;
            int top_;

            constraint_board_t
            create_constraint_board(const std::vector<const constraint_t *> &constraints) {
                
                // the constraint board cb holds a collection cb[i] of to-be-checked constraints
                // for each variable i.
                constraint_board_t cb(vars_.size());

                // for each constraint
                for( const auto &c : constraints ) {

                    // find 'last' variable of constraint c
                    // by iterating over the variables of the constraint 
                    int last_idx = -1;
                    for ( auto var : c->vars() ) {
                        // skip determined variables of the constraint
                        if ( a_[var] != Undetermined ) continue;

                        // determine maximum index of this variable (in the order of the variables
                        // of the assignment iterator)
                        //    determine the index of the constraint variable in the iterator variables
                        auto idx = std::distance(vars_.begin(), find(vars_.begin(),vars_.end(),var));
                        
                        //    determine maximum as last_idx
                        last_idx = std::max(last_idx, static_cast<int>(idx));
                    }

                    if ( 0 <= last_idx && last_idx < static_cast<int>(vars_.size()) ) {
                        cb[last_idx].push_back( c );
                    }
                    // std::cout << "Constraint on ";
                    // std::copy(c->vars().begin(),c->vars().end(),std::ostream_iterator<int>(std::cout,","));
                    // std::cout<<" last idx: "<<last_idx<<" of ";
                    // std::copy(vars_.begin(),vars_.end(),std::ostream_iterator<int>(std::cout,", "));
                    // for (auto var : c->vars()) { 
                    //     if ( a_[var] != Undetermined ) 
                    //         std::cout << "det " << var <<" = " << a_[var] << " ";
                    // }
                    // std::cout << std::endl;

                }

                return cb;
            }
        }; // end AssignmentIterator

        template<class Constraint>
        auto
        make_iterator(const std::vector<var_idx_t> &vars,
                      const std::vector<int> &domsizes,
                      const std::vector<const Constraint *> &constraints) {
            return AssignmentIterator<Constraint>(*this, vars, domsizes, constraints);
        }

    private:
        std::vector<var_value_t> values_;
    };

    const int Assignment::Undetermined;
}

#endif
