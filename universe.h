#pragma once 

#include "randomfield.h"
#include "background.h"
#include "main.h"
#include "pcg32.h"
#include "timer.h"

#include <set>
#include <vector>
#include <algorithm>
#include <memory>
#include <string>


class Universe { 

    friend class Multiverse;

public:

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
    CollisionTasks::iterator leftBoundaryTask;


    /* the factor of 1/8: it so happens that 4 particles per fourier bin of the initial spectrum avoids seeing top-hat related artifacts. but we don't want such high-frequency IC components in the system
     * beyond the ones that can be resolved by the particles - it will be aliasing noise at much lower frequencies otherwise. Factor 1/2 may be enough. But if we want to restrict to 4 particles / fourier bin max fourier frequency for FFT analysis, we should have another 1/4... */
    Universe(Background bg, Float L, const Long seed, Long nParticles, int boundaryCondition) : L(L), pcg(0, seed), nParticles(nParticles), initDisplacement(nParticles/2), bg(bg), Dx(L/(nParticles-1)), boundary(boundaryCondition) {
        //std::cout << "Bang! ";
        if (!bg.isIntegrated) {
            std::cout << "Integrating background (because it has not yet been done)." << std::endl; 
            this->bg.integrate(); 
        }
    }

    // force move construction, not only for performance (e.g. because of particles std::vector):
    // we cannot copy the collisions std::multiset tree and expect the iterator on it in Particle to still be valid
    // this crashes the evolution after Multiverse::bang() ! 
    Universe& operator=(const Universe&) = delete;
    Universe(const Universe&) = delete;
    Universe& operator=(Universe&&) = default;
    Universe(Universe&&) = default;

    ~Universe() {};

    void draw(PowerSpectrum* ps) {
        initDisplacement.generate(pcg, ps);
        sampleParticles();
    }

    Float get_particle_pos_standardized(Long i) const  { return (particles[i].x + Dx/2 + L/2)/(L+Dx); }
    Float get_real_particle_pos_standardized(Long i) const  { return (particles[i].x + Dx/2 + L/2)/(L+Dx) + particles[i].sheet; }
    //Float get_particle_pos(Long i) const  { return particles[i].x + particles[i].sheet*L; }
    Float get_particle_vel(Long i) const  { return particles[i].v; }
    Float get_particle_time(Long i) const { return particles[i].t; }
    Long get_particle_id(Long i) const { return particles[i].index; }
    Long get_particle_collision_number(Long i) const { return particles[i].nCollisions; }
    ParticleList get_particles() const {return particles;}

    const Long nParticles; 
    const Float L; 
    const Float Dx;

    static constexpr int PERIODIC = 1;
    static constexpr int REFLECTIVE = 2;
    int boundary;
    int rotation = 0;

    void update_collision();

    void synchronize() {
        for (int i = 0; i < nParticles; ++i) 
            update_particle(i, most_recent_particle_time());
    }

    void evolve(Float z, bool skipCollisions = false) {

        Float t = bg.getTauOfZ(z);

        if (!skipCollisions) { 
            while (next_collision_time() < t)
                update_collision();
        }

        latestTime = t;
        synchronize();
    }

    void diagnose(std::string m = "") {
        //for (Long i = 0; i < nParticles-1; ++i) { 
            //if (particles[i+1].x < particles[i].x && std::fabs(particles[i+1].t - particles[i].t) < 1e-15){
        //std::cout << "Diagnose with L = " << L << std::endl;
        //std::cout << m << std::endl;
                //std::cout << "Particle (" << i << ", " <<  i + 1 << ") are in the wrong order." << std::endl;
                //std::cout << "Coords are " << particles[i].x << ", " << particles[i+1].x << std::endl;
                //assert(false);
            //}
            //if (particles[i].x < -L/2-Dx/2  || particles[i].x > L/2 + Dx/2 ) { 
        //std::cout << "Diagnose with L = " << L << std::endl;
        //std::cout << m << std::endl;
                //std::cout << "Particle " << i << " is out of the box at " << particles[i].x << std::endl;
                //assert(false);
            //}
        //}
    }

    // returns time of last collision or current evolved time if newer 
    Float most_recent_particle_time() const { return latestTime; }

    Float next_collision_time() const {return (collisions.begin())->collTime; }

private:

    long nCollisions = 0;
    // update particle to time t assuming no collision
    void update_particle(Long i, Float t);


    void update_particle_task(Long i) {
        //Float collTime = particles[i].t + collision_time(i);
        //std::cout << "attempting computation of coll time for " << i << std::endl;
        
        //if (i==9) {
            //if (particles[9].t > 20.343139812720) {
                //std::cout << "Found moment " << std::endl;
                //std::cout << "particle is at " << particles[9].x << " " << particles[9].v << std::endl;
                ////update_particle(9,  25.83775386);
                //std::cout << "particle will be at " << particles[9].x << std::endl;
            //}
        //}
        //
        //START_NAMED_TIMER("COLLTIME")
        Float collTime = collision_time(i);
        //STOP_TIMER
        //std::cout << "Coll time just computed for " << i <<" is " << collTime << std::endl;
        //
        //START_NAMED_TIMER("INSERT")
        if (i >= 0) { 
            particles[i].task = collisions.insert(RightCollision(i-rotation, collTime));
        
            
            if (collTime <= particles[i].t) 
                std::cout << "Bad time order at " << i << " with " << (collTime-particles[i].t)/collTime << std::endl; 
            
            assert(collTime > particles[i].t);
        }
        else 
            leftBoundaryTask = collisions.insert(RightCollision(i-rotation, collTime));
        //STOP_TIMER
    }

    Float latestTime = 0;
    // compute right-collision time of partile idLeft
    Float collision_time(Long idLeft);
    //SomeArray<Particle> 
    // particle list must be sorted by position
    //


    pcg32 pcg;

    void sampleParticles(); 

public:
    Background bg; 
    RandomField initDisplacement;
private:
    Float find_max_timestep();
    void integrate_step(Float timestep);
};

