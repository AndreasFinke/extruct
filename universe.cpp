#include "universe.h"

#include <iostream>

void Universe::sampleParticles() {

    particles.clear();
    particles.reserve(nParticles);

    // fill particles according to displacement field 
    Float totalVel = 0;
    for (Long i = 0; i < nParticles; ++i) {
        Float x = (Float(i)+Float(0.5)) / nParticles;
        //std::cout << x << std::endl;
        Float disp = initDisplacement.get_field(x);
        //std::cout << disp << " "; 
        // something like \dot(H) needed here
        Float dh = 1;
        Float vel  = disp*dh;
        particles.push_back(Particle(x+disp, vel));
        //
        totalVel += vel;
        //initial time
        particles[i].t = 0;
    }

    // set total momentum to zero
    totalVel /= nParticles;
    for (Long i = 0; i < nParticles; ++i) 
        particles[i].v -= totalVel;

    // ensure order
    bool unsorted = false;
    for (Long i = 1; i < nParticles; ++i) {
        if (particles[i].x < particles[i-1].x) {
            unsorted = true;
            break;
        }
    }
    if (unsorted)  {
        std::cout << "Initial state contains stream crossing. Sorting particlces..." << std::endl;
        std::sort(particles.begin(), particles.end(), [](const Particle& lhs, const Particle& rhs) { return lhs.x < rhs.x;} );
        std::cout << "... done." << std::endl;
    }

    // tasks 
    collisions.clear();
    for (Long i = 0; i < nParticles-1; ++i) 
        update_particle_task(i);

    particles[nParticles-1].task = collisions.end();

}

//Float force(Long id, Long n, Float x) { 
Float force(Long id, Long n) { 
    return Float(0.01)*((n-1)*Float(0.5)-id);
    //return Float(0.01)*((n-1)*Float(0.5)-id + x-Float(0.5));
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
    //if (vd > 0)
        return vd + sqrt(vd*vd+2*xd*F); //TODO improve numerics
    //else
        //return (-2*xd*F)/(vd - sqrt(vd*vd+2*xd*F));
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

    Float t = coll->collTime;
    latestTime = t;
    Long id = coll->id;
    
    // std::cout << "coll part is " << id << std::endl;

    // this collision is not needed anymore. 

    collisions.erase(coll);

    // update two colliding particles
    // (the _tasks_ of particle at id (still pointing to coll) and neighbors will be updated below)

    update_particle(id, t);
    update_particle(id+1, t);

    // positions should now agree - check 
    assert(std::fabs(particles[id].x - particles[id+1].x) < 0.00000001);

    // swap colliding particle positions in sorted particle list - that's in their future 
   
    std::swap(particles[id], particles[id+1]);


    // if we were not at boundary and (id+2) exists then update it 
    // We know there cannot have been another collision by time t!
    if (id + 1 < nParticles-1) {
        update_particle(id + 2, t);
        // all up to 3 particles 
        //    (id) <- collided -> (id+1)   (id+2)
        // involved are now at the same time t, 
        // so we can use collision_time without problem 
        // note also that all other collision times further out 
        // (e.g. (id+2) (id+3) ) are unaffected by what happened here (i.e. id+2 moving to time t)   
        // their collision times with their neighbors were absolute times and forces do not change
        
        //update the collision to (id+2)
        update_particle_task(id+1); 
        // there was a nontrivial collision with id+2 previously, which is now at particle id (due to swap) - delete it! 
        collisions.erase(particles[id].task);
    }
    else { // case that id+1 is the boundary
        // set the new end
        particles[id+1].task = collisions.end();
        // the old no-task end() is at particle[id].task and needs no delete
        // (will be overwritten below)
    }
    // if id-1 exists, update its collision with new id (which was id+1) 
    if (id > 0) {
        // bring id-1 to time t as well
        update_particle(id-1, t);
        // old collision invalid, remove
        collisions.erase(particles[id-1].task);
        // id-1 and id are at time t, can use: 
        update_particle_task(id-1);
    }

    // finally, and in any case, update necessarily nontrivial task for particle id 
    update_particle_task(id);

}
