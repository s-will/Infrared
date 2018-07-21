#include <iostream>
#include <vector>
#include "infrared.hpp"
#include "rnadesign.hpp"

using namespace ired;

int
main(int argc, char **argv) {
    Assignment a(3);
    std::vector<int> vars {0,1,2};
    std::vector<int> domsizes {4,4,4};
    auto cc1 = ComplConstraint(0,1);
    auto cc2 = ComplConstraint(0,2);
    auto constraints = std::vector<const Constraint*> {&cc1,&cc2};

    std::vector<int> vars0 {};
    auto constraints0 = std::vector<const Constraint*> {};

    for( auto it0 = a.make_iterator(vars0, domsizes, constraints0); ! it0.finished(); ++it0 ) {
        for( auto it = a.make_iterator(vars, domsizes, constraints); ! it.finished(); ++it ) {
            std::copy(a.values().begin(),a.values().end(),std::ostream_iterator<int>(std::cout," "));
            std::cout << std::endl;
        }
    }
}
