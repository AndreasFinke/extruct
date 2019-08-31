#pragma once

#include "universe.h"
#include "background.h"
#include "main.h"
#include "measurement.h"

#include "tbb/tbb.h"

#include <vector>
#include <memory>


class Multiverse {

public:

    Multiverse() {}

    int bang(Long nParticles, const Background& bg, PowerSpectrum* pspec, Float L) {
        universes.push_back(Universe(bg, pspec, L, universes.size(), nParticles));
        return universes.size()-1;

    }

    int nUniverses() { return universes.size(); }

    
    void evolve(int universeID, Float z) {

        universes[universeID].evolve(z);

    }

    void evolveAll(Float z) {
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, universes.size(), 1),
                [&](tbb::blocked_range<size_t> range) {
                    for (size_t i = range.begin(); i < range.end(); ++i) {
                    //for (int i = 0; i < universes.size(); ++i) {
                        std::cout << "evolving universe " << i << " until z = " << z << std::endl;
                        universes[i].evolve(z);
                        nCollisions += universes[i].nCollisions;
                        std::cout << "... done after " << universes[i].nCollisions << " collisions." << std::endl; 
                        universes[i].nCollisions = 0;
                    }
                    std::cout << "Total number of collisions in this multiverse: " << nCollisions << std::endl; 
                }
        );

    }

    void measure(int universeID, Measurement* m) {
        assert(universeID < universes.size());
        m->measure(universes[universeID], 1);
    }

    void measureAll(Measurement* m) {
        
        for (auto& u : universes) 
            m->measure(u, universes.size());
    }

private:

    std::vector<Universe> universes;
    long nCollisions = 0;

};
