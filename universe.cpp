#include "universe.h"

#include <iostream>

void Universe::sampleParticles() {

    particles.clear();
    particles.reserve(nParticles);

    // fill particles according to displacement field 
    Float totalVel = 0;
    for (Long i = 0; i < nParticles; ++i) {
        Float x = i*Dx;
        //std::cout << x << std::endl;
        /* random field expects 0..1 grid coords */
        Float disp = initDisplacement.get_field(x/L);
        x -= L/2;
        //std::cout << disp << " "; 
        // something like \dot(H) needed here
        Float dh = 1;
        Float vel  = disp*dh;

        if (x+disp < -L/2) 
            x += L;
        else if (x+disp > L/2)
            x -= L;

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
    for (Long i = -1; i < nParticles; ++i) 
        update_particle_task(i);

    //particles[nParticles-1].task = collisions.end();

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
//Float Universe::collision_time(Long idLeft) {
    //assert( std::fabs(particles[idLeft].t - particles[idLeft+1].t) < 0.00000000001 );
    //Float vd = particles[idLeft+1].v - particles[idLeft].v;
    //Float xd = particles[idLeft+1].x - particles[idLeft].x;
    //Float F = force(1, 2) - force(2, 2);
    //F = 1/F;
    //vd *= F;
    ////if (vd > 0)
        //return vd + sqrt(vd*vd+2*xd*F); //TODO improve numerics
    ////else
        ////return (-2*xd*F)/(vd - sqrt(vd*vd+2*xd*F));
//}


/*returns the absolute collision time of particle idLeft with its right neighbor. For -1, computes considers collision of first particle with left boundary at -L/2. For nParticles-1, considers collision of last particle with boundary at L/2 */

Float Universe::collision_time(Long idLeft) {
    Float Delta0, Deltad0, tau0;
   
    if (idLeft == -1) {
        Delta0 = particles[0].x + L/2;
        Deltad0 = particles[0].v;
        tau0 = particles[0].t;
    }
    else if (idLeft == nParticles-1) {
        Delta0 = L/2 - particles[idLeft].x;
        Deltad0 = - particles[idLeft].v;
        tau0 = particles[idLeft].t;
    }
    else {
        Delta0  = particles[idLeft+1].x - particles[idLeft].x;
        Deltad0 = particles[idLeft+1].v - particles[idLeft].v;
        tau0 = particles[idLeft].t; 
        assert( std::fabs(tau0 - particles[idLeft+1].t) < 0.00000000001 );
    }

    assert( Delta0 > -0.0000001 );

    //Float a0 = bg.getScaleFactor(tau0);
    Float D10 = bg.getD1(tau0);
    Float D20 = bg.getD2(tau0);
    Float D1d0 = bg.getD1d(tau0);
    Float D2d0 = bg.getD2d(tau0);
    Float Pec0 = bg.getPec(tau0)*Dx;
    Float Pecd0 = bg.getPecd(tau0)*Dx;

    Float c1 = - (Delta0 - Pec0)*D2d0 + (Deltad0 - Pecd0)*D20;
    Float c2 =   (Delta0 - Pec0)*D1d0 - (Deltad0 - Pecd0)*D10;

    int idxLast = bg.NTABLE - 1;

    auto distance = [c1, c2, Dx=Dx, &bg=bg] (int idx) { return c1*bg.D1[idx] + c2*bg.D2[idx] + bg.Pec[idx]*Dx; } ;
    auto distanced = [c1, c2, Dx=Dx, &bg=bg] (int idx) { return c1*bg.D1d[idx] + c2*bg.D2d[idx] + bg.Pecd[idx]*Dx; } ;

    /* there is at most one zero crossing into the negative in the future. if the function is positive at the end, there is no zero or a zero is in the future after the final time considered. 
       we just return the final time plus one plus the particle index here */
    if (distance(idxLast) > 0)
        return bg.taufin + 1 + idLeft;
    
    /* else, find sign change in the future */

    int idxLeft = std::min(int(tau0/bg.dtau), idxLast-1);

    /* add one to be in the future of the start time */
    int idxR = idxLeft + 1;
    /* walk until sign changed; note it will happen at the latest at idxLast because of the previous if statement*/
    for (; distance(idxR) > 0; ++idxR);
    /* the zero is contained in the interval formed with the previous index */ 
    int idxL = idxR - 1;

    Float yL = distance(idxL);
    Float yR = distance(idxR);
    Float DyL = distanced(idxL);
    Float DyR = distanced(idxR);

    Float mL = bg.getOm()*1.5*bg.a[idxL];
    Float fL = -Dx*mL;
    Float mR = bg.getOm()*1.5*bg.a[idxR];
    Float fR = -Dx*mR;
    Float DDyL = mL*yL + fL;
    Float DDyR = mR*yR + fR;

    return Spline::find_zero(yL, DyL, DDyL, yR, DyR, DDyR, idxL*bg.dtau, bg.dtau);

}


void Universe::update_particle(Long id, Float tau) {
    Float fac = (id - (nParticles-1)*Float(0.5))*Dx;
    Float tau0 = particles[id].t;
    //Float T = tau - tau0;
    particles[id].t =  tau;

    Float a0 = bg.getScaleFactor(tau0);
    Float D10 = bg.getD1(tau0);
    Float D20 = bg.getD2(tau0);
    Float D1d0 = bg.getD1d(tau0);
    Float D2d0 = bg.getD2d(tau0);
    Float Pec0 = bg.getPec(tau0)*fac;
    Float Pecd0 = bg.getPecd(tau0)*fac;

    Float c1 = - (particles[id].x - Pec0)*D2d0 + (particles[id].v - Pecd0)*D20;
    Float c2 =   (particles[id].x - Pec0)*D1d0 - (particles[id].v - Pecd0)*D10;

    int idxLast = bg.NTABLE - 1;

    //std::cout << "Particle " << id << " was at " << particles[id].x << " and is now (tau= " << tau << ") at ";
    particles[id].x = c1*bg.getD1(tau) + c2*bg.getD2(tau) + bg.getPec(tau)*fac; 
    //std::cout << particles[id].x << std::endl;
    particles[id].v = c1*bg.getD1d(tau) + c2*bg.getD2d(tau) + bg.getPecd(tau)*fac; 

    //particles[id].x += particles[id].v * T + Float(0.5)*T*T*force(id, nParticles);
    //particles[id].v += T * force(id, nParticles); 
}

void Universe::update_collision() {

    // pick next collision

    auto coll = collisions.begin();
    // now coll == particles[coll->id].task unless it is the left boundary collision 

    Float t = coll->collTime;
    latestTime = t;
    Long id = coll->id;
    
     //std::cout << "coll part is " << id << std::endl;
    
    nCollisions++; 

    // this collision is not needed anymore. 

    collisions.erase(coll);

    // update two colliding particles
    // (the _tasks_ of particle at id (still pointing to coll) and neighbors will be updated below)

    if (id == -1) {
        update_particle(0, t);
        //std::cout << "After boundary collision first particle is at " << particles[0].x << std::endl;
        particles[0].v = -particles[0].v;
        update_particle_task(-1);
        collisions.erase(particles[0].task);
        update_particle(1, t);
        update_particle_task(0);
    }
    else if (id == nParticles-1) { 
        update_particle(id, t);
        //std::cout << "After boundary collision last particle is at " << particles[id].x << std::endl;
        particles[id].v = -particles[id].v;
        update_particle_task(id);
        collisions.erase(particles[id-1].task);
        update_particle(id-1, t);
        update_particle_task(id-1);
    }
    else {

        update_particle(id, t);
        update_particle(id+1, t);

        // positions should now agree - check 
        //std::cout << "pos should now agree - id " << id << " error " << std::fabs(particles[id].x - particles[id+1].x)/particles[id].x << std::endl; 
        //assert(std::fabs(particles[id].x - particles[id+1].x) < 0.00000001);

        // swap colliding particle positions in sorted particle list - that's in their future 
      
        Float meanpos =  (particles[id].x + particles[id+1].x)*Float(0.5);
        particles[id].x = particles[id+1].x = meanpos;
            
        std::swap(particles[id], particles[id+1]);


        // if we were not at boundary and (id+2) exists then update it 
        // We know there cannot have been another collision by time t!
        if (id + 1 < nParticles-1) {
            update_particle(id + 2, t);
            if (particles[id+2].x <= meanpos) {
                std::cout.precision(18);
                std::cout << "oops. " << id + 2 << " is at " << particles[id+2].x << " and meanpos is " << meanpos << std::endl;
                particles[id+2].x = meanpos + 0.00000001;
                std::cout << "now " << id + 2 << " is at " << particles[id+2].x << std::endl;
            }
            assert(particles[id+2].x > meanpos);
            // the particles particles 
            //    (id) <- collided and swapped -> (id+1)   (id+2)
            // are now at the same time t, 
            // so we can use collision_time called by update_particle_task without problem and find new right collision time of (id+1) (id+2)
            // note also that all other collision times further out 
            // (e.g. (id+2) (id+3) ) are unaffected by what happened here (i.e. id+2 moving to time t)   
            // their collision times with their neighbors were absolute times and forces do not change
            
            //update the collision of (id+1) with (id+2)
            update_particle_task(id+1); 
            // there was a nontrivial collision with id+2 previously, which is now at particle id (due to swap) - delete it! 
            collisions.erase(particles[id].task);
        }
        else { // case that id+1 is the last particle
             //set the new no-task end() to the boundary particle 
            //particles[id+1].task = collisions.end();
            /* compute new collision with right boundary */
            update_particle_task(id+1);
            // from before, the task with boundary collision is assigned to particle[id].task, but needs no delete
            // (will be overwritten below)
        }
        // if id-1 exists, update its collision with new id (which was id+1) 
        if (id > 0) {
            // bring id-1 to time t as well
            update_particle(id-1, t);
            if (particles[id-1].x >= meanpos) {
                std::cout.precision(18);
                std::cout << "oops. " << id - 1 << " is at " << particles[id-1].x << " and meanpos is " << meanpos << std::endl;
                particles[id-1].x = meanpos - 0.00000001;
                std::cout << "now " << id - 1 << " is at " << particles[id-1].x << std::endl;
            }
            assert(particles[id-1].x < meanpos);
            // old collision invalid, remove
            collisions.erase(particles[id-1].task);
            // id-1 and id are at time t, can use: 
            update_particle_task(id-1);
        }
        else {
            update_particle_task(-1);
        }

        // finally, and in any case, update necessarily nontrivial task for particle id 
        update_particle_task(id);

    }

}
