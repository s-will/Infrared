#ifndef INFRARED_RNADESIGN_HPP
#define INFRARED_RNADESIGN_HPP

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
 * @brief RNA design-specific extension of Infrared
 *
 * This file provides domain-specific InfraRed extensions for RNA design sampling
 */


#include "iostream"

namespace ired {

    namespace rnadesign {
        
        /**
         * @brief Base pair complementarity constraint
         *
         * Constrains two variables to represent complementary base
         * pairs. Due to A=0,C=1,G=2,U=3, this limits the value
         * combinations to AU=03, CG=12, GU=23 and symetrically UA=30,
         * GC=21, UG=32.
         */
        class ComplConstraint : public Constraint {
        public:
            using self_t = ComplConstraint;
            using parent_t = Constraint;
            using base_t = typename parent_t::base_t;
            
            /**
             * @brief Construt to constrain two variables to be
             * complementary
             *
             * @param i index of first variable
             * @param j index of second variable
             */
            ComplConstraint(int i, int j)
                : Constraint({i,j}) {
            }
            
            /**
             * @brief Evaluate constraint on {i,j}
             *
             * @param a assignment
	     * @return bool
             */
            bool
            operator ()(const Assignment &a) const override {

                static std::array<bool,16>
                    tab = {false,false,false,true , //A
                           false,false,true ,false, //C
                           false,true ,false,true , //G
                           true ,false,true ,false};//U

                auto is_compl = tab[ a[vars()[0]] + 4 * a[vars()[1]] ];
                // std::cout << "Check "
                //           << vars()[0] << " " 
                //           << vars()[1] << " " 
                //           << a[vars()[0]] << " " 
                //           << a[vars()[1]] << " " 
                //           << is_compl << std::endl;
                return
                    is_compl;
            }
        };

        /**
         * @brief Same complementarity class constraint
         *
         * Ensure that bases are in the same compoment of the bipartition
         * due to complementarity
         * Components are {AG} and {CU}
         */
        class SameComplClassConstraint : public Constraint {
        public:
            using self_t = ComplConstraint;
            using parent_t = Constraint;
            using base_t = typename parent_t::base_t;

            SameComplClassConstraint(int i, int j)
                : Constraint({i,j}) {
            }

            bool
            operator ()(const Assignment &a) const override {
                auto x = a[vars()[0]];
                auto y = a[vars()[1]];
                return (x&1) == (y&1);
            }
        };
    
        /**
         * @brief Different complementarity class constraint
         *
         * Ensure that bases are *not* in the same compoment of the
         * bipartition due to complementarity. @see SameComplClassConstraint
         */
        class DifferentComplClassConstraint : public Constraint {
        public:
            using parent_t = Constraint;
            using base_t = typename parent_t::base_t;

            DifferentComplClassConstraint(int i, int j)
                : Constraint({i,j}) {
            }

            bool
            operator ()(const Assignment &a) const override {
                auto x = a[vars()[0]];
                auto y = a[vars()[1]];
                return (x&1) != (y&1);
            }
        };

        /**
         * @brief Function for control of GC content
         */
        class GCControl : public Function<double> {
        public:
            using parent_t = Function<double>;
            using base_t = typename parent_t::base_t;

	    /**
	     * @brief Constructor of GC content control
	     * @param i index of variable
	     * @param weight GC ratio targeting value
	     */
            GCControl(int i , double weight)
                : parent_t({i}), weight_(weight) {
            }

	    /**
	     * @brief Evaluate GC content on position i
	     * @return weight if assignment on position is C or G,
	     * otherwise 1.
	     */
            double
            operator ()(const Assignment &a) const override {

                static std::array<double,4> tab = {1,weight_,weight_,1};

                return
                    tab[ a[vars()[0]] ];
            }

        private:
            double weight_;
        };
  
        /**
         * @brief Base class of RNA energy functions
         *
         * Supports static polymorphism (CRTP) to simplify the
         * definition of RNA energy functions.
         */
        template<class T>
        class RNAEnergyFunction : public Function<double> {
        public:
            using parent_t = Function<double>;

            RNAEnergyFunction( std::vector<typename parent_t::var_idx_t> vars, double weight)
                : parent_t(vars), weight_(weight) {
            }

	    /**
	     * @brief Energy evaluation of the assignment
	     * @param a assignment
	     * @return 
             * @f[ \pi^{E(a)} @f]
	     */
            double
            operator ()(const Assignment &a) const override {
                return
                    pow ( weight_, - static_cast<const T*>(this)->energy(a) );
            }
        
        protected:

            //! index of bp: in order AU, UA, CG, GC, GU, UG;
            //! index=index of outer + 6 * index of inner
            static
            int
            bpindex(int x, int y) {
                return bpindex_tab_[x+4*y];
            }

            double weight_;

            //! table storing index of bp
            static std::array<int, 16> bpindex_tab_;
        };

        /**
         * @brief Base pair energy
         *
         * Distinguishs stacked and terminal bps.
         * Holds static table of base pair energy parameters.
         */
        class BPEnergy : public RNAEnergyFunction<BPEnergy> {
        public:
            using parent_t = RNAEnergyFunction;
            using base_t = typename parent_t::base_t;

            BPEnergy(int i, int j, double weight)
                : parent_t({i,j} , weight), is_terminal_(false) {
            }
            
            BPEnergy(int i, int j, bool is_terminal, double weight)
                : parent_t({i,j}, weight), is_terminal_(is_terminal) {
            }

            /** @brief set energy table
             *
	     * @param table a vector of size 6
             * order of the 6 parameters in the table: AU,CG,GU,AU_term,CG_term,GU_term
             */
            static
            void
            set_energy_table(const std::vector<double> table) {
                assert(table.size() == 6);
                tab_.resize(6);
                std::copy(table.begin(),table.end(),tab_.begin());
            }

            //! @brief get energy
            double 
            energy(const Assignment &a) const {
                if (tab_.size()!=6) {
                    std::cerr<<"Table of stacking energies has to be initialized before use."<<std::endl;
                    exit(-1);
                }

                auto i = vars()[0];
                auto j = vars()[1];
                auto x = bpindex(a[i],a[j]);
            
                if (x<0) {
                    return std::numeric_limits<double>::infinity();
                }
                x=x/2;
            
                return
                    tab_[ x + (is_terminal_?3:0) ];
            }
            
        private:
            static std::vector<double> tab_;
            bool is_terminal_;
        };

        /**
         * @brief Stacking energy
         *
         * Holds static table of stacking energy parameters.
         */
        class StackEnergy : public RNAEnergyFunction<StackEnergy> {
        public:
            using parent_t = RNAEnergyFunction;
            using base_t = typename parent_t::base_t;

            //! @brief construct energy function for stack
            //! of base pairs (i,j) with (i+1,j-1)
            StackEnergy(int i, int j, double weight)
                : parent_t({i,j,i+1,j-1}, weight) {
            }

            //! @brief set energy table for stacks
            //!
            static
            void
            set_energy_table(const std::vector<double> table) {
                assert(table.size()==18);
                tab_.resize(18);
                std::copy(table.begin(),table.end(),tab_.begin());
            }

            //! @brief get energy
            double
            energy(const Assignment &a) const {
                if (tab_.size()!=18) {
                    std::cerr<<"Table of stacking energies has to be initialized before use."<<std::endl;
                    exit(-1);
                }

                auto i = vars()[0];
                auto j = vars()[1];

                auto x = bpindex(a[i],a[j]);
                auto y = bpindex( a[i+1], a[j-1] );
            
                if (x<0 || y<0) {
                    return std::numeric_limits<double>::infinity();
                }
            
                if (x & 1) {
                    std::swap(x,y);
                    y ^= 1; //xor -> swap base pair order (AU->UA, etc)
                }

                return tab_[ y + 6*(x/2) ];
            }

        private:
            static std::vector<double> tab_;
        };

        // define static members
        
        std::vector<double> BPEnergy::tab_;
        std::vector<double> StackEnergy::tab_;

        //! index 0..5 in order AU, UA, CG, GC, GU, UG
        template<class T>
        std::array<int,16> RNAEnergyFunction<T>::
        bpindex_tab_ = { -1,-1,-1, 0,  //A
                         -1,-1, 2,-1,  //C
                         -1, 3,-1, 4,  //G
                          1,-1, 5,-1 };//U

    } // end namespace rnadesign
} //end namespace ired
#endif

