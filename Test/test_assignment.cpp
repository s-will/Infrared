#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#undef NDEBUG

#include <iostream>
#include <vector>

#include "../Infrared/infrared.hpp"
#include "../Infrared/rnadesign.hpp"
#include "../Infrared/functions.hpp"

using namespace ired;
using namespace ired::rnadesign;

using EvaluationPolicy = StdEvaluationPolicy<double>;

// wrap a constraint as double valued function

class ConstraintFunction : public Function<double> {
public:
    ConstraintFunction(const Constraint *c) : Function<double>(c->vars()), c_(c) {}

    double operator() (const Assignment &a) const override {
        return (*c_)(a)?1.0:0.0;
    }

private:
    const Constraint *c_;
};

int
main(int argc, char **argv) {
    BPEnergy::set_energy_table({-2,-3,-1,-2,-3,-1});

    // Test standard case of enumeration
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

    int counter=0;

    std::cout <<"Enum some assignments; consistent constraints" << std::endl;

    for( auto it = a.make_iterator(vars, cn, constraints, functions, 1.0);
         ! it.finished(); ++it ) {
        std::copy(a.values().begin(),a.values().end(),
                  std::ostream_iterator<int>(std::cout," "));
        std::cout <<":" << it.value() << std::endl;
        counter++;
    }
    assert(counter == 14);

    // Test enumeration over inconsistent constraints
    auto cc4 = ComplConstraint(1,3);
    constraints.push_back( &cc4 );
    counter=0;

    std::cout <<"Try enumeration; inconsistent constraints" << std::endl;

    for( auto it = a.make_iterator(vars, cn, constraints, functions, 1.0);
         ! it.finished(); ++it ) {
        std::copy(a.values().begin(),a.values().end(),
                  std::ostream_iterator<int>(std::cout," "));
        std::cout <<":" << it.value() << std::endl;
        counter++;
    }
    assert(counter == 0 );

    std::cout <<"Try enumeration; inconsistency with functions" << std::endl;

    constraints.pop_back();
    auto f3 = ConstraintFunction(&cc4);
    functions.push_back( &f3 );

    for( auto it = a.make_iterator(vars, cn, constraints, functions, 1.0);
         ! it.finished(); ++it ) {
        std::copy(a.values().begin(),a.values().end(),
                  std::ostream_iterator<int>(std::cout," "));
        std::cout <<":" << it.value() << std::endl;
        counter++;
    }
    assert(counter == 0 );

    // Test enumeration over empty variables
    counter=0;
    std::vector<int> vars0 {};
    auto constraints0 = std::vector<const Constraint*> {};
    auto functions0 = std::vector<const Function<double>*> {};

    for( auto it = a.make_iterator(vars0, cn, constraints0, functions0, 1.0);
         ! it.finished(); ++it ) {
        counter++;
    }

    assert(counter == 1);
}
