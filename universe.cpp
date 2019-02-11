#pragma once 

#include "universe.h"

#include <iostream>

void Universe::sampleParticles() {

    particles.clear();
    particles.reserve(nParticles);

    // fill particles according to displacement field 
    for (Long i = 0; i < nParticles; ++i) {
        Float x = (Float(i)+Float(0.5)) / nParticles;
        //std::cout << x << std::endl;
        Float disp = initDisplacement.get_field(x)/nParticles;
        // something like \dot(H) needed here
        Float dh = 1;
        Float vel  = disp*dh;
        particles.push_back(Particle(x+disp, vel));

        //initial time
        particles[i].t = 0;
    }

    // confirm order 
    bool flag = false;
    for (Long i = 1; i < nParticles; ++i) {
        if (particles[i].x < particles[i-1].x) {
            flag = true;
            break;
        }
    }
    if (flag) 
        std::cout << "WARNING: initial state contains stream crossing. Particles are not sorted." << std::endl;

    // tasks 
    collisions.clear();
    for (Long i = 0; i < nParticles-1; ++i) {
        Float collTime = particles[i].t + collision_time(i); 

        // create new collision task for this particle and store pointer to the task in the particle 
        particles[i].task = collisions.insert(RightCollision(i, collTime));
    }
    particles[nParticles-1].task = collisions.end();


}

Float force(Long id, Long n) { 
    return Float(0.01)*(n*Float(0.5)-id);
}

/* compute right-collision times of particle at idLeft with particle at idLeft+1
 * assuming both particles are at the same time
 * assuming constant force \ddot{x} = Force(x) = nParticles/2 - index(x)
 * xL + vL T + F(xL) T^2 / 2 = xR + vR T  + F(xR) T^2/2 
 * F T^2 + 2(vL-vR) T - 2(xR-xL) = 0
 * T = (vR-vL)/F +- sqrt((vR-vL)^2/F^2 + 2(xR-xL)/F). negative sign is collision in the past, positive in the future
 *
 * note: the rightmost particle is ignored and task list will have one element less than number of particles. 
 * It may still collide and its left neighbor may become the rightmost particle. 
*/
Float Universe::collision_time(Long idLeft) {
    assert( std::fabs(particles[idLeft].t - particles[idLeft+1].t) < 0.00000000001 );
    Float vd = particles[idLeft+1].v - particles[idLeft].v;
    Float xd = particles[idLeft+1].x - particles[idLeft].x;
    Float F = force(1, 2) - force(2, 2);
    F = 1/F;
    vd *= F;
    return vd + sqrt(vd*vd+2*xd*F); //TODO improve numerics
}


void Universe::update_particle(Long id, Float t) {
    Float T = t - particles[id].t;
    particles[id].t =  t;

    particles[id].x += particles[id].v * T + Float(0.5)*T*T*force(id, nParticles);
    particles[id].v += T * force(id, nParticles); 
}

void Universe::update_collision() {
   
    // pick next collision

    auto coll = collisions.begin();
    assert(coll == particles[coll->id].task);
    
    std::cout << "coll part is " << coll->id << std::endl;
    // update two colliding particles
   
    Float t = coll->collTime;
    latestTime = t;
    Long id = coll->id;

    update_particle(id, t);
    update_particle(id+1, t);

    //positions should now agree - check 

    assert(std::fabs(particles[id].x - particles[id+1].x) < 0.00000001);
    std::cout << std::fabs(particles[coll->id].x - particles[coll->id+1].x) << std::endl;

    // swap their positions in sorted particle list since that's their future 
   
    std::swap(particles[id], particles[id+1]);

    // update the (up to) two surrounding particles. We know they cannot have collided by time t!
    if (id + 1 < nParticles-1) {
        update_particle(id + 2, t);
    }

    // all up to 3 particles 
    //    (id) <- collided -> (id+1)   (id+2)
    // involved are now at the same time t, 
    // so we can use collision_time without problem 
    // note also that all other collision times further out 
    // (e.g. (id+2) (id+3) ) are unaffected by what happened here (i.e. id+2 moving to time t)   
    // their collision times with their neighbors were absolute times and forces do not change
    
    auto update_task = [&](Long id) { 
        collisions.erase( particles[id].task ); 
        Float collTime = particles[id].t + collision_time(id);
        particles[id].task = collisions.insert(RightCollision(id, collTime));
    };

    update_task(id); // note: coll has been deleted now!

    if (id + 1 < nParticles-1)
        update_task(id+1);
    else {
        particles[id+1].task = collisions.end();
    }
}
