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
 */


namespace ired {
    class ComplConstraint : public Constraint {
    public:
        using self_t = ComplConstraint;
        using parent_t = Constraint;
        using base_t = typename parent_t::base_t;

        ComplConstraint(int i, int j)
            : Constraint({i,j}) {
        }

        ComplConstraint(const self_t &c)
            : Constraint(c) {
        }

        bool
        operator ()(const Assignment &a) const override {

            static std::array<bool,16>
                tab = {0,0,0,1, //A
                       0,0,1,0, //C
                       0,1,0,1, //G
                       1,0,1,0};//U

            return
                tab[ a[vars()[0]] + 4 * a[vars()[1]] ];
        }

        std::unique_ptr<base_t>
        clone() const override {
            return std::make_unique<self_t>(*this);
        }
    };

    class CGControl : public Function<double> {
    public:
        using self_t = CGControl;
        using parent_t = Function<double>;
        using base_t = typename parent_t::base_t;

        CGControl(int i , double *weight)
            : parent_t({i}), weight_(weight) {
        }

        CGControl(const self_t &c)
            : parent_t(c), weight_(c.weight_) {
        }

        double
        operator ()(const Assignment &a) const override {

            static std::array<double,4> tab = {0,-1,-1,0};

            return
                pow ( *weight_, - tab[ a[vars()[0]] ] );
        }

        std::unique_ptr<base_t>
        clone() const override {
            return std::make_unique<self_t>(*this);
        }

    private:
        double *weight_;
    };


    class BPEnergy : public Function<double> {
    public:
        using self_t = BPEnergy;
        using parent_t = Function<double>;
        using base_t = typename parent_t::base_t;

        BPEnergy(int i, int j, double *weight)
            : parent_t({i,j}), weight_(weight) {
        }

        BPEnergy(const self_t &c)
            : parent_t(c), weight_(c.weight_) {
        }

        double
        operator ()(const Assignment &a) const override {

            static std::array<double,16> tab =
                {0,0,0,-2, //A
                 0,0,-3,0, //C
                 0,-3,0,-1, //G
                 -2,0,-1,0};//U

            return
                pow ( *weight_, - tab[ a[vars()[0]] + 4 * a[vars()[1]] ] );
        }

        std::unique_ptr<base_t>
        clone() const override {
            return std::make_unique<self_t>(*this);
        }

    private:
        double *weight_;
    };
}

#endif
