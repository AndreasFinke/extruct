#pragma once 

#include "randomfield.h"
#include "background.h"
#include "main.h"
#include "pcg32.h"

#include <set>
#include <vector>
#include <algorithm>
#include <memory>


class Universe { 

    friend class Multiverse;

public:

    Universe(const Background& bg, PowerSpectrum* ps, Float L, const Long seed, Long nParticles) : L(L), pcg(0, seed), initDisplacement(512), nParticles(nParticles), bg(bg), Dx(L/(nParticles-1)) {
        std::cout << "Universe ctor. " << bg.isIntegrated << std::endl; 
        draw(ps);
    }

    // force move construction, not only for performance (e.g. because of particles std::vector):
    // we cannot copy the collisions std::multiset tree and expect the iterator on it in Particle to still be valid
    // this crashes the evolution after Multiverse::bang() ! 
    Universe& operator=(const Universe&) = delete;
    Universe(const Universe&) = delete;
    Universe(Universe&&) = default;

    ~Universe() {};

    void draw(PowerSpectrum* ps) {
        initDisplacement.generate(pcg, ps);
        sampleParticles();
    }

    Float get_particle_pos_standardized(Long i) const  { return particles[i].x/L + 0.5; }
    Float get_particle_pos(Long i) const  { return particles[i].x; }
    Float get_particle_vel(Long i) const  { return particles[i].v; }
    Float get_particle_time(Long i) const { return particles[i].t; }

    const Long nParticles; 
    const Float L; 
    const Float Dx;


    void update_collision();

    void synchronize() {
        for (int i = 0; i < nParticles; ++i) 
            update_particle(i, most_recent_particle_time());
    }

    void evolve(Float t) {
        while (next_collision_time() < t)
            update_collision();

        latestTime = t;
        synchronize();
    }

    // returns time of last collision or current evolved time if newer 
    Float most_recent_particle_time() { return latestTime; }

    Float next_collision_time() {return (collisions.begin())->collTime; }

private:

    long nCollisions = 0;
    // update particle to time t assuming no collision
    void update_particle(Long i, Float t);


    void update_particle_task(Long i) {
        //Float collTime = particles[i].t + collision_time(i);
        Float collTime = collision_time(i);
        particles[i].task = collisions.insert(RightCollision(i, collTime));
    }

    Float latestTime = 0;
    // compute right-collision time of partile idLeft
    Float collision_time(Long idLeft);
    //SomeArray<Particle> 
    // particle list must be sorted by position
    //

    Background bg; 

    pcg32 pcg;


    struct RightCollision {
        //RightCollision(RightCollision&&) = default;
        //RightCollision(const RightCollision&) = default;
        RightCollision(Long particleId, Float collisionTime) : id(particleId), collTime(collisionTime) {} 
        /* the idx-th particle in  to collide with its right neighbor */  
        Long id = 0;
        /* ... in time T */
        Float collTime = 0;
        bool operator<(const RightCollision& rhs) const {return collTime < rhs.collTime;}
    };

    // the right datastructure for a priority queue is a heap. think of a binary tree with the only condition that both children of any node have a value that is larger than the node. it's only partially sorted and fast to retrieve the top element and insert new ones. 
    // however we must be able to also modifiy some elements (spatial neighbor particles) as we go so the heap is unfortunately out:
    // std has a prioriry_queue heap that does not provide functionality beyond retrieving the top, and a heapified vector but it would be a mess to deal with that manually (and perhaps not be guaranteed to run on all std implementations) ((and suffer from the same problems as the sorted vector below))
    // the next good way for a priority queue is a fully sorted data structure 
    // think of a sorted array, where we can find in logarithmic time and access the top in constant time
    // finding in log time needs random access, so lists are out (they are anyway slow in practice). a sorted vector works, but to insert an element a swap is not enough - lots of things have to be rotated by one position! I tried and it becomes slow if the inserts are not often in the early <10000 elements which seems to have to do with cache size. I don't think we can guarantee that. One might think deques would be good because they should consist of shorter chunks of contiguous memory but in practice they are even slower
    // a set / multiset is internally sorted. it works very well and beats vector unless inserts are early. Then the sorted vector is up to 10-20 times faster. but for late inserts the set stayes very fast and beats the sorted vector by any factor 
    // on top of that, as mentioned, we must access the neighbor particles etc.. were we to use a vector, we would keep indices into the task list. but once a particle is re-inserted into the sorted list, not only do we need to copy around part of the sorted list, we also must update the corresponding indices of the particles into the task list and there seems no good way to do so - one would have to keep pointers to the particles in each task list ! to be able to update the index in the particle... so find the host particle, find its neighbors, find their tasks, see if they have an updated index, okay. but also iterate through ALL tasks with rotated index and find all corresponding particles and update their tasks... a huge mess, and potentially much slower 
    // with a set, instead of indices into the set we keep iterators (pointers!). each particle gets an iterator on its task list set element. these are preserved when a new element is inserted into the set!!! 
    // internally, a set is some tree (typically red-black) which is a lot of unexposed functionality guaranteeing it works well! neat, that we can use that here. it is precisely the pointer-to-node functionality that vector does not have and the fast insert that list cannot provide that is both available here. 
    //
    using CollisionTasks = std::multiset<RightCollision>; 
    CollisionTasks collisions; 
  
    using Particle = TaskParticle<CollisionTasks::iterator>;
    using ParticleList = std::vector<Particle>;

    ParticleList particles; 

    void sampleParticles(); 

    RandomField initDisplacement;
    Float find_max_timestep();
    void integrate_step(Float timestep);
};

