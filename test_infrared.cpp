///////////////////////////////////////////////////////////
//
// Test usage of sample.so in C++
//
////////////////////////////////////////////////////////////

#include <memory>
#include <iostream>
#include <array>
#include <ctime>

#include "infrared.hpp"
#include "rnadesign.hpp"

using namespace ired;


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

    srand(time(0)); // just kidding ;)


    int seqlen=30;

    auto ctd = ClusterTree<>(seqlen, 4);

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
