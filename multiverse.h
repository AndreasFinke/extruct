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

    Multiverse(int seed = 0, bool onlyZeldovich = false) : seed(seed), onlyZeldovich(onlyZeldovich) {}

    int bang(Long nParticles, const Background& bg, PowerSpectrum* pspec, Float L, int boundaryCondition) {
        std::cout << "Drawing universe " << universes.size() << " with seed " << universes.size()+seed << "\n";
        universes.push_back(Universe(bg, L, universes.size() + seed, nParticles, boundaryCondition));
        universes.back().draw(pspec);
        return universes.size()-1;

    }

    void bangAll(Long nParticles, const Background& bg, PowerSpectrum* pspec, Float L, int boundaryCondition, int N) {
        for (int i = 0; i < N; ++i) {
            universes.push_back(Universe(bg, L, universes.size() + seed, nParticles, boundaryCondition));
        }
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, universes.size(), 1),
                [&](tbb::blocked_range<size_t> range) {
                    for (size_t i = range.begin(); i < range.end(); ++i) {
                    //for (int i = 0; i < universes.size(); ++i) {
                    
                        universes[i].draw(pspec);
                        std::cout << "Drawing universe " << i << " with seed " << i+seed << "\n";
                    }
                }
        );

    }

    int nUniverses() { return universes.size(); }

    
    void evolve(int universeID, Float z) {

        universes[universeID].diagnose();
        universes[universeID].evolve(z, onlyZeldovich);

    }

    void evolveAll(Float z) {
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, universes.size(), 1),
                [&](tbb::blocked_range<size_t> range) {
                    for (size_t i = range.begin(); i < range.end(); ++i) {
                    //for (int i = 0; i < universes.size(); ++i) {
                        std::cout << "evolving universe " << i << " until z = " << z << std::endl;
                        evolve(i, z);
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
        
        for (auto&& u : universes) 
            m->measure(u, universes.size());
        m->postprocess();
    }

private:

    std::vector<Universe> universes;
    Long nCollisions = 0;
    Long seed = 0;

    bool onlyZeldovich; 

};
