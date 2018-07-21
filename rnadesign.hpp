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
 *
 * This file provides domain-specific InfraRed extensions for RNA design sampling
 *
 * @todo find a python compatible way to allow later changes of weights without redefining everything
 * ( needs indirection; similar to passing a pointer to weight instead of passing it directly )
 */


#include "iostream"

namespace ired {
    class ComplConstraint : public Constraint {
    public:
        using self_t = ComplConstraint;
        using parent_t = Constraint;
        using base_t = typename parent_t::base_t;

        ComplConstraint(int i, int j)
            : Constraint({i,j}) {
        }

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

    class GCControl : public Function<double> {
    public:
        using self_t = GCControl;
        using parent_t = Function<double>;
        using base_t = typename parent_t::base_t;

        GCControl(int i , double weight)
            : parent_t({i}), weight_(weight) {
        }

        double
        operator ()(const Assignment &a) const override {

            static std::array<double,4> tab = {0,-1,-1,0};

            return
                pow ( weight_, - tab[ a[vars()[0]] ] );
        }

    private:
        double weight_;
    };


    class BPEnergy : public Function<double> {
    public:
        using self_t = BPEnergy;
        using parent_t = Function<double>;
        using base_t = typename parent_t::base_t;

        BPEnergy(int i, int j, double weight)
            : parent_t({i,j}), weight_(weight) {
        }

        static
        void
        set_energy_table(const std::vector<double> table) {
            assert(table.size() == 16);
            std::copy(table.begin(),table.end(),tab_.begin());
        }

        double
        operator ()(const Assignment &a) const override {
            return
                pow ( weight_, - tab_[ a[vars()[0]] + 4 * a[vars()[1]] ] );
        }

    private:
        double weight_;

        static std::array<double,16> tab_;
    };
    
   #define INF 1.0e6
   std::array<double,16> BPEnergy::
   tab_ = { INF,INF,INF, -2, //A
            INF,INF, -3,INF, //C
            INF, -3,INF, -1, //G
             -2,INF, -1,INF};//U
    #undef INF
   
}



#endif
