#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#undef NDEBUG

#include <iostream>
#include <vector>

#include "../Infrared/infrared.hpp"
#include "../Infrared/rnadesign.hpp"

using namespace ired;
using namespace ired::rnadesign;

using EvaluationPolicy = StdEvaluationPolicy<double>;

int
main(int argc, char **argv) {
    BPEnergy::set_energy_table({-2,-3,-1,-2,-3,-1});

    Assignment a(4);
    std::vector<int> vars {0,1,2,3};
    std::vector<int> domsizes {4,4,4,4};
    
    auto cn = ConstraintNetwork<double>(domsizes);
    
    auto cc1 = ComplConstraint(0,1);
    auto cc2 = ComplConstraint(1,2);
    auto cc3 = ComplConstraint(0,3);
    auto f1 = GCControl(0, 0.1);
    auto f2 = BPEnergy(2, 3, 2);
    auto constraints = std::vector<const Constraint*> {&cc1,&cc2,&cc3};
    auto functions = std::vector<const Function<double>*> { &f1, &f2 };

    std::vector<int> vars0 {};
    auto constraints0 = std::vector<const Constraint*> {};
    auto functions0 = std::vector<const Function<double>*> {};
    
    int counter=0;
    
    for( auto it = a.make_iterator(vars, cn, constraints, functions, 1.0); 
         ! it.finished(); ++it ) {
        std::copy(a.values().begin(),a.values().end(),
                  std::ostream_iterator<int>(std::cout," "));
        std::cout <<":" << it.value() << std::endl;
        counter++;
    }

    assert(counter == 14);
    
    counter=0;
    
    for( auto it = a.make_iterator(vars0, cn, constraints0, functions0, 1.0); 
         ! it.finished(); ++it ) {
        counter++;
    }
    
    assert(counter == 1);
}
