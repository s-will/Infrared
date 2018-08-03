///////////////////////////////////////////////////////////
//
// Test usage of infrared in C++
//
////////////////////////////////////////////////////////////

#include <memory>
#include <iostream>
#include <array>
#include <ctime>

#include "../Infrared/infrared.hpp"
#include "../Infrared/rnadesign.hpp"

using namespace ired;
using namespace ired::rnadesign;

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

    srand(time(0)); // poor man's seeding

    // initialize static energy table for all BPEnergy functions
    BPEnergy::set_energy_table({-2,-3,-1,-2,-3,-1});

    int seqlen=30;
    int num_seqs=20;

    auto ct = ClusterTree<>(seqlen, 4);

    // ((((.)())).(.))
    auto basepairs = std::vector<std::pair<int,int>>{ {1,9}, {2,8}, {3,5}, {6,7}, {11,13}, {0,14} };

    auto added = std::vector<bool>(seqlen);

    double weightE = 5.0;
    double weightE2 = 10.0;
    double weightGC = 0.5;

    // add clusters for base pairs
    for (auto p : basepairs) {
        // build Cluster

        auto cluster = ct.add_root_cluster( {p.first,p.second} );
        ct.add_constraint( cluster, std::make_unique<ComplConstraint>(p.first, p.second) );
        ct.add_function( cluster, std::make_unique<BPEnergy>(p.first, p.second, weightE ) );
        ct.add_function( cluster, std::make_unique<GCControl>(p.first, weightGC ) );
        ct.add_function( cluster, std::make_unique<GCControl>(p.second, weightGC ) );

        added[p.first] = true;
        added[p.second] = true;

        // hack in some extra constraints! (second target structure)
        if (p.second==14) {
            auto cluster2 = ct.add_child_cluster( cluster, {14,20} );
            ct.add_function( cluster2, std::make_unique<GCControl>(20, weightGC ) );
            ct.add_constraint( cluster2, std::make_unique<ComplConstraint>(14,20) );
            ct.add_function( cluster2, std::make_unique<BPEnergy>(14, 20, weightE2) );
            added[20] = true;

            auto cluster3 = ct.add_child_cluster( cluster, {14,25,27,29} );
            // 14 25
            ct.add_function( cluster3, std::make_unique<GCControl>(25, weightGC ) );
            ct.add_constraint( cluster3, std::make_unique<ComplConstraint>(14,25) );
            ct.add_function( cluster3, std::make_unique<BPEnergy>(14, 25, weightE2) );
            added[25] = true;

            // 25 27
            ct.add_function( cluster3, std::make_unique<GCControl>(27, weightGC ) );
            ct.add_constraint( cluster3, std::make_unique<ComplConstraint>(25,27) );
            ct.add_function( cluster3, std::make_unique<BPEnergy>(25, 27, weightE2) );
            added[27] = true;

            // 27 29
            ct.add_function( cluster3, std::make_unique<GCControl>(29, weightGC ) );
            ct.add_constraint( cluster3, std::make_unique<ComplConstraint>(29,27) );
            ct.add_function( cluster3, std::make_unique<BPEnergy>(29, 27, weightE2) );
            added[29] = true;

            // 14 29
            ct.add_function( cluster3, std::make_unique<GCControl>(29, weightGC ) );
            ct.add_constraint( cluster3, std::make_unique<ComplConstraint>(14,29) );
            ct.add_function( cluster3, std::make_unique<BPEnergy>(14, 29, weightE2) );

            // generate a further child cluster with more overlap --> message size 4^3
            auto cluster4 = ct.add_child_cluster( cluster3, {14,25,28,29} );
            added[28] = true;
            ct.add_function( cluster4, std::make_unique<GCControl>(28, weightGC ) );
            ct.add_constraint( cluster4, std::make_unique<ComplConstraint>(25,28) );
            ct.add_function( cluster4, std::make_unique<BPEnergy>(25, 28, weightE2) );
            ct.add_constraint( cluster4, std::make_unique<ComplConstraint>(28,29) );
            ct.add_function( cluster4, std::make_unique<BPEnergy>(28, 29, weightE2) );
        }

        if (p.second==13) {
            auto cluster2 = ct.add_child_cluster( cluster, {13,21} );
            ct.add_function( cluster2, std::make_unique<GCControl>(21, weightGC ) );
            ct.add_constraint( cluster2, std::make_unique<ComplConstraint>(13,21) );
            ct.add_function( cluster2, std::make_unique<BPEnergy>(13, 21, weightE2) );
            added[21] = true;
        }

    }

    // add clusters singletons
    for (int i=0; i<seqlen; i++) {
        if (added[i]) continue;

        auto cluster = ct.add_root_cluster({i});
        ct.add_function( cluster, std::make_unique<GCControl>(i, weightGC ) );
    }

    // write dot bracket
    std::vector<char> s(seqlen,'.');
    for (auto p : basepairs) {
        s[p.first] = '(';
        s[p.second] = ')';
    }
    std::copy(s.begin(),s.end(),std::ostream_iterator<char>(std::cout,""));
    std::cout<<std::endl;

    ct.evaluate();

    for (int i=0; i<num_seqs; i++) {
        const Assignment &a = ct.sample();
        std::cout << a << std::endl;
    }
}
