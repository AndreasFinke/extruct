#pragma once

#include "universe.h"
#include "background.h"
#include "main.h"
#include "measurement.h"


#include <vector>
#include <memory>


class Multiverse {

public:

    Multiverse() {}

    int bang(Long nParticles, const Background& bg, PowerSpectrum* pspec) {
        universes.push_back(Universe(bg, pspec, universes.size(), nParticles));
        return universes.size()-1;

    }

    int nUniverses() { return universes.size(); }

    
    void evolve(int universeID, Float t) {

        universes[universeID].evolve(t);

    }

    void evolveAll(Float t) {

        for (int i = 0; i < universes.size(); ++i) {
            std::cout << "evolving universe " << i << " until t=" << t << std::endl;
            universes[i].evolve(t);
        }

    }

    void measure(int universeID, Measurement* m) {
        m->measure(universes[universeID], 1);
    }

    void measureAll(Measurement* m) {
        
        for (auto u : universes) 
            m->measure(u, universes.size());
    }

private:

    std::vector<Universe> universes;

};
