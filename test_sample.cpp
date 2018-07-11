///////////////////////////////////////////////////////////
//
// Test usage of sample.so in C++
//
////////////////////////////////////////////////////////////

#include "sample.hpp"
#include <memory>
#include <iostream>
#include <array>

using namespace bscf;

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



std::ostream &
operator <<(std::ostream &out, const Assignment &a) {
    auto num2nucl = [](int x){ return std::string("ACGU")[x]; };
    auto vals = a.values();
    std::string seq;
    std::transform(vals.begin(),vals.end(),std::back_inserter(seq),num2nucl);
    return out << seq;
}




int
main(int argc, char** argv) {
    int seqlen=30;

    auto ctd = ClusterTreeDecomposition<>(seqlen, 4);

    // ((((.)())).(.))
    auto basepairs = std::vector<std::pair<int,int>>{ {1,9}, {2,8}, {3,5}, {6,7}, {11,13}, {0,14} };

    auto added = std::vector<bool>(seqlen);

    double weightE = 10.0;
    double weightCG = 0.5;

    // add clusters for base pairs
    for (auto p : basepairs) {
        // build Cluster

        auto cluster = ctd.make_root_cluster( {p.first,p.second} );
        ctd.add_constraint( cluster, ComplConstraint(p.first, p.second) );
        ctd.add_function( cluster, BPEnergy(p.first, p.second, &weightE ) );
        ctd.add_function( cluster, CGControl(p.first, &weightCG ) );
        ctd.add_function( cluster, CGControl(p.second, &weightCG ) );
        
        added[p.first] = true;
        added[p.second] = true;
    }

    // add clusters singletons
    for (int i=0; i<seqlen; i++) {
        if (added[i]) continue;

        auto cluster = ctd.make_root_cluster({i});
        ctd.add_function( cluster, CGControl(i, &weightCG ) );
    }

    // write dot bracket
    std::vector<char> s(seqlen,'.');
    for (auto p : basepairs) {
        s[p.first] = '('; 
        s[p.second] = ')';
    }
    std::copy(s.begin(),s.end(),std::ostream_iterator<char>(std::cout,""));
    std::cout<<std::endl;

    ctd.evaluate();

    for (int i=0; i<10; i++) {
        const Assignment &a = ctd.sample();
        std::cout << a << std::endl;
    }
}
