#include "universe.h"

#include <iostream>

void Universe::sampleParticles() {

    particles.clear();
    particles.reserve(nParticles);

    // fill particles according to displacement field 
    Float totalVel = 0;
    for (Long i = 0; i < nParticles; ++i) {

        Float x = (i+1)*L/(nParticles+1);
        x = i*L/(nParticles-1);


        /* random field expects 0..1 grid coords and does not include a length dimension*/
        Float disp = L*initDisplacement.get_displacement(x/L);
        x -= L/2;


        /* if displ = D1*A, vel is D1d*A = D1d/D1*displ */
        /* but for now my D1 is 0 initially... the growing mode of matter would instead be the scale factor. 
         * Then we just find the Hubble rate (which is D2) in place of the ratio! -  But with a factor of a^2 for getting the tau-derivative! */

        if (bg.zin < 10 && bg.Om < 0.9) 
            std::cout << "Initial redshift under 10 in a non - Einstein - de Sitter - univere. Initial particle veloctities may be inaccurate." << std::endl;
        if (std::fabs(bg.Oc)/(bg.zin+1) > 0.001 ) 
            std::cout << "Nonzero curvature at late starting redshift causes an error of about " << 100*std::fabs(bg.Oc)/(bg.zin+1) << " \% on the initial particle velocities." << std::endl;

        Float vel =  disp*bg.h*bg.D2[0]*bg.ain*bg.ain;

        short sheet = 0;
        if (x+disp < -L/2 - Dx/2)  {
            std::cout << "Warning: Large initial displacments caused wrap-around at boundary." << std::endl; 
            if (boundary == REFLECTIVE) {
                disp += 2*(-L/2-Dx/2-x-disp);
                vel = -vel;
            }
            else {
                x += L+Dx;
                sheet--;
            }
        }
        else if (x+disp > L/2 + Dx/2) { 
            std::cout << "Warning: Large initial displacments caused wrap-around at boundary." << std::endl; 
            if (boundary == REFLECTIVE) {
                disp -= 2*(x+disp-L/2-Dx/2);
                vel = -vel;
            }
            else { 
                x -= (L+Dx);
                sheet++;
            }
        }

        particles.push_back(Particle(x+disp, vel, i, sheet));
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

            //for (int i = 0; i < nParticles; ++i) std::cout << " " << particles[i].index;
            //std::cout << std::endl << std::endl;
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
    Float extraDelta = 0; 
    Float fac = Dx;
    Float order = 1;   
   
    if (idLeft == -1) {
        //fac *= 0.5; //TODO
        fac *= (0.5-rotation+nParticles*(particles[idLeft+1].sheet));
        extraDelta = particles[0].sheet*(L+Dx);
        Delta0 = particles[0].x + extraDelta - (-L/2-Dx/2); //TODO
        Deltad0 = particles[0].v;
        tau0 = particles[0].t;
    //assert( Delta0 >= -1e-15 );
    }
    else if (idLeft == nParticles-1) {
        //fac *= 0.5; //TODO
        fac *= (0.5+rotation+nParticles*( - particles[idLeft].sheet));
        extraDelta = -particles[idLeft].sheet*(L+Dx);
        Delta0 = (L/2+Dx/2) - particles[idLeft].x + extraDelta; //TODO?
        Deltad0 = -particles[idLeft].v;
        tau0 = particles[idLeft].t;
        //std::cout << particles[idLeft].x << " and sheet " << particles[idLeft].sheet << " extraD " << extraDelta << " Delta0 " << Delta0 << std::endl;
    //assert( Delta0 >= -1e-15 );
    }
    else {

        fac *= (1 + nParticles*(particles[idLeft+1].sheet - particles[idLeft].sheet));//TODO
        //std::cout << particles[idLeft].x << " " << particles[idLeft+1].x;
        
        extraDelta = (particles[idLeft+1].sheet - particles[idLeft].sheet)*(L+Dx);
        Delta0  = particles[idLeft+1].x - particles[idLeft].x + extraDelta; //TODO
        //Delta0  = particles[idLeft+1].x - particles[idLeft].x;
        Deltad0 = particles[idLeft+1].v - particles[idLeft].v;
        tau0 = particles[idLeft].t; 

        //std::cout << " at " << tau0 << std::endl;
        assert( std::fabs(tau0 - particles[idLeft+1].t) < 1e-15 );
    }

    if (Delta0 - extraDelta <= -1e-15) {
        if (idLeft >= -1)
            std::cout <<"Left x " << particles[idLeft].x << " and sheet " << particles[idLeft].sheet << std::endl;
        if (idLeft + 1 < nParticles)
            std::cout << "Right x " << particles[idLeft+1].x << " and sheet " << particles[idLeft +1].sheet <<  std::endl;
        std::cout << "Delta0 " << Delta0 << " extraDelta " << extraDelta << "Diff " << Delta0 - extraDelta << std::endl;
        std::cout << "idLeft " << idLeft << " sheet L " << particles[idLeft >= 0 ? idLeft : idLeft + 1].sheet << std::endl;
    }
    assert(Delta0 - extraDelta > -1e-15);

    //if (std::fabs(extraDelta) > L/2 ) { 
        //std::cout << "Potential ordering issue, sign of extraDelta is " << sgn(extraDelta)  << std::endl;
        //order = sgn(Delta0 + 1e-15); 
        //Delta0 *= order;
        //Deltad0 *= order;
        //fac *= order;
        //std::cout << "order is " << order << std::endl;
    //}
    

    //std::cout << Delta0 << std::endl;

    //assert( Delta0 >= -1e-15 );

    //Float a0 = bg.getScaleFactor(tau0);
    Float D10 = bg.getD1(tau0);
    Float D20 = bg.getD2(tau0);
    Float D1d0 = bg.getD1d(tau0);
    Float D2d0 = bg.getD2d(tau0);
    Float Pec0 = bg.getPec(tau0)*fac;
    Float Pecd0 = bg.getPecd(tau0)*fac;

    Float c1 = - (Delta0 - Pec0)*D2d0 + (Deltad0 - Pecd0)*D20;
    Float c2 =   (Delta0 - Pec0)*D1d0 - (Deltad0 - Pecd0)*D10;

    int idxLast = bg.NTABLE - 1;

    auto distance = [c1, c2, fac, extraDelta, &bg=bg] (int idx) { return c1*bg.D1[idx] + c2*bg.D2[idx] + bg.Pec[idx]*fac - extraDelta; } ; //TODO 
    auto distanced = [c1, c2, fac, &bg=bg] (int idx) { return c1*bg.D1d[idx] + c2*bg.D2d[idx] + bg.Pecd[idx]*fac; } ;

    /* there is at most one zero crossing into the negative in the future in the standard case (fac > 0). if the function is positive at the end, there is no zero or a zero is in the future after the final time considered. 
       we just return the final time plus one plus the particle index here */

    /* else, find sign change in the future */

    int idxLeft = std::min(int(tau0/bg.dtau), idxLast-1);

    /* add one to be in the future of the start time */
    int idxR = idxLeft + 1;

    if (fac >= 0) {
        if (distance(idxLast) > 0)
            return bg.taufin + 2 + idLeft;
        /* walk until sign changed; note it will happen at the latest at idxLast because of the previous if statement*/
        for (; distance(idxR) > 0; ++idxR);
    }
    else {
        if (Deltad0 > 0)
            return bg.taufin + 2 + idLeft;
        else {
            for (; distance(idxR) > 0; ++idxR) {
                if (idxR == idxLast)
                    return bg.taufin + 2 + idLeft;
            }
        }
    }
    
    /* the zero is contained in the interval formed with the previous index */ 
    int idxL = idxR - 1;

    Float yL = distance(idxL);
    Float yR = distance(idxR);
    Float DyL = distanced(idxL);
    Float DyR = distanced(idxR);

    Float mL = bg.getOm()*1.5*bg.a[idxL];
    Float fL = -fac*mL;   /// TODO 
    Float mR = bg.getOm()*1.5*bg.a[idxR];
    Float fR = -fac*mR;
    Float DDyL = mL*yL + fL;
    Float DDyR = mR*yR + fR;

    if (!std::isnormal(DDyL) || !std::isnormal(DDyR) ) {
        std::cout << "Anormal second derivatives in collision_time " << idxL << " of " << idxLast <<  std::endl;
        std::cout << "tau0 " << tau0 << " idxLeft " << idxLeft << " distance(idxLeft+1) " << distance(idxLeft + 1) << std::endl;
        std::cout << "c1 " << c1 << " c2 " << c2 << " bg.D1 " << bg.D1[idxLeft+1] << " bg.D2d " << bg.D2d[idxLeft+1] << std::endl;
        std::cout << "D10 " << D10 << " D20 " << D20 << " D1d0 " << D1d0 << " D2d0 " << D2d0 << " Pec0 " << Pec0 << " Pecd0 " << Pecd0 << " fac " << fac << " Delta0 " << Delta0 << " Deltad0 " << Deltad0 << std::endl;
        std::cout << "idLeft " << idLeft << " x " << particles[idLeft].x << " x2 " << particles[idLeft+1].x << " extraD " << extraDelta << std::endl;
        std::cout << bg.getOm() << " " <<  bg.a[idxL] << std::endl;
    }

    return Spline::find_zero(yL, DyL, DDyL, yR, DyR, DDyR, idxL*bg.dtau, bg.dtau);

}

void Universe::update_particle(Long id, Float tau) {
    Float fac = (id - rotation + particles[id].sheet*nParticles - (nParticles-1)*Float(0.5))*Dx;
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

    Float c1 = - (particles[id].x + particles[id].sheet * (L + Dx) - Pec0)*D2d0 + (particles[id].v - Pecd0)*D20;
    Float c2 =   (particles[id].x + particles[id].sheet * (L + Dx) - Pec0)*D1d0 - (particles[id].v - Pecd0)*D10;

    int idxLast = bg.NTABLE - 1;

    //std::cout << "Particle " << id << " was at " << particles[id].x << " and is now (tau= " << tau << ") at ";
    particles[id].x = c1*bg.getD1(tau) + c2*bg.getD2(tau) + bg.getPec(tau)*fac - particles[id].sheet * (L + Dx); 
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
    Long id = coll->id + rotation;
   
    //Long id = coll->id;

    //for (int i = 0; i < nParticles; ++i) std::cout << " " << particles[i].index<< " " << particles[i].t;
    //std::cout << std::endl<<"coll part is " << id << " at " << t;
    //if (id >= 0) 
        //std::cout << " with index " << particles[id].index << std::endl;
    //else
        //std::cout << std::endl;
    
    nCollisions++; 

    // this collision is not needed anymore. 

    collisions.erase(coll);

    // update two colliding particles
    // (the _tasks_ of particle at id (still pointing to coll) and neighbors will be updated below)

    std::cout.precision(18);
    if (id == -1) {
        collisions.erase(particles[0].task);
        update_particle(0, t);
        //std::cout << "L After boundary collision first particle is at " << particles[0].x << std::endl;
        //
        if (boundary == REFLECTIVE) { 
            particles[0].x = -L/2 - Dx/2;
            particles[0].v = -particles[0].v;
            update_particle_task(-1);
            update_particle(1, t);
            update_particle_task(0);
        }
        else if (boundary == PERIODIC) {

            /* beam particle to the right boundary*/
            particles[0].x = L/2 + Dx/2;
            
             //std::cout << "Now at  " << particles[0].x << std::endl;
            std::cout << "L";
            /* remember where it came from */ 
            --particles[0].sheet;

            //for (auto particle : particles) std::cout << particle.index << " ";
            //std::cout << std::endl;
            /* put this previously leftmost particle at the right end of memory in sorted list as well */
            std::rotate(particles.begin(), std::next(particles.begin()), particles.end());
            //for (auto particle : particles) std::cout << particle.index << " ";
            //std::cout << std::endl;

            rotation--;

            //for (CollisionTasks::iterator it = collisions.begin(); it != collisions.end(); ++it)
                //const_cast<Long&>((*it).id)--;

            //for (CollisionTasks::iterator it = collisions.begin(); it != collisions.end(); ++it)
                //std::cout << (*it).id << " ";
            //std::cout << std::endl; 
            
            //for (int i = 0; i < nParticles; ++i) std::cout << " " << particles[i].index<< "(" << particles[i].sheet << ") "  << particles[i].t;
            //std::cout << std::endl;

            /* update left boundary collision with new leftmost particle */
            update_particle_task(-1);

            /* advance previous rightmost particle to current time and forget its old collision scheduled */
            update_particle(nParticles-2, t);
            collisions.erase(particles[nParticles-2].task);

            /* find new collision time of previous rightmost and current rightmost particle */
            update_particle_task(nParticles-2);

            /* find new collision time of current rightmost particle and right boundary */
            update_particle_task(nParticles-1);
            //for (CollisionTasks::iterator it = collisions.begin(); it != collisions.end(); ++it)
                //std::cout << (*it).id << " ";
            //std::cout << std::endl; 
        }

        diagnose("left boundary");
    }
    else if (id == nParticles-1) { 
        collisions.erase(particles[id-1].task);
        //std::cout << "Before boundary collision last particle is at " << particles[id].x << std::endl;
        update_particle(id, t);
        //std::cout << "R After boundary collision last particle is at " << particles[id].x << std::endl;
        //std::cout << "Coll time was " << t << std::endl;
        std::cout << "R"; 
        //
        if (boundary == REFLECTIVE) { 
            particles[id].x = L/2 + Dx/2;
            particles[id].v = -particles[id].v;
            update_particle_task(id);
            update_particle(id-1, t);
            update_particle_task(id-1);
        }
        else if (boundary == PERIODIC) {

            /* beam to left boundary */
            particles[nParticles-1].x = -L/2 - Dx/2;

            /* remember origin */
            ++particles[nParticles-1].sheet;

            /* take previous rightmost as new beginning of the array (new leftmost) */
            //for (auto particle : particles) std::cout << particle.index << " ";
            //std::cout << std::endl;
            std::rotate(particles.begin(), std::prev(particles.end()), particles.end());
            //for (auto particle : particles) std::cout << particle.index << " ";
            //std::cout << std::endl;
            rotation++;
            //for (int i = 0; i < nParticles; ++i) std::cout << " " << particles[i].index<< "(" << particles[i].sheet << ") "  << particles[i].t;
            //std::cout << std::endl;
            //for (CollisionTasks::iterator it = collisions.begin(); it != collisions.end(); ++it)
                //const_cast<Long&>((*it).id)++;

            //for (CollisionTasks::iterator it = collisions.begin(); it != collisions.end(); ++it)
                //std::cout << (*it).id << " ";
            //std::cout << std::endl; 
            //for (auto co : collisions) {
                //co.id++;
            //}

            /* update right boundary collision */
            update_particle_task(nParticles-1);
            /* take previous leftmost to current time */
            update_particle(1, t);
            /* and find collision with it */
            update_particle_task(0);
            /* update left boundary collision */
            collisions.erase(leftBoundaryTask);
            update_particle_task(-1);
            //for (CollisionTasks::iterator it = collisions.begin(); it != collisions.end(); ++it)
                //std::cout << (*it).id << " ";
            //std::cout << std::endl; 
        }

        diagnose("right boundary");
    }
    else {

        //std::cout << "Index of collision particle is " << particles[id].index << std::endl;
        update_particle(id, t);
        update_particle(id+1, t);
        particles[id].nCollisions++;
        particles[id+1].nCollisions++;

        // positions should now agree - check 
        //std::cout << "pos should now agree - id " << id << " error " << std::fabs(particles[id].x - particles[id+1].x)/particles[id].x << std::endl; 
        if (std::fabs(particles[id].x - particles[id+1].x) > 1e-8*L)
        {
            std::cout.precision(18);
            std::cout << "Unequal positions of " << id << ", " << id+1 << " at times " << particles[id].t << " " << particles[id+1].t <<  std::endl;
            std::cout << "These carry index " << particles[id].index << ", " << particles[id+1].index << " and sheets " << particles[id].sheet << " " << particles[id+1].sheet << std::endl;
            std::cout << "Coords are " << particles[id].x << ", " << particles[id+1].x << " and rotation is " << rotation << std::endl;
            assert(false);
        }

        // swap colliding particle positions in sorted particle list - that's in their future 
        std::swap(particles[id], particles[id+1]);

        Float meanpos =  (particles[id].x + particles[id+1].x)*Float(0.5);

        particles[id].x = meanpos;
        particles[id+1].x = meanpos;
       // particles[id].x -= 1e-15;
       // particles[id+1].x += 1e-15;
        diagnose("after swapping");


        /* if id+2 is a real particle */ 
        if (id + 1 < nParticles-1) {
            update_particle(id + 2, t);

            if (particles[id+2].x <= meanpos) {
                std::cout.precision(18);
                std::cout << "oops. " << id + 2 << " is at " << particles[id+2].x << " and meanpos is " << meanpos << std::endl;
                particles[id+2].x = meanpos + 0.00000001;
                std::cout << "now " << id + 2 << " is at " << particles[id+2].x << std::endl;
            }
            assert(particles[id+2].x > meanpos);
        }
        // the particles particles 
        //    (id) <- collided and swapped -> (id+1)   (id+2)
        // are now at the same time t, 
        // so we can use collision_time called by update_particle_task without problem and find new right collision time of (id+1) (id+2)
        // note also that all other collision times further out 
        // (e.g. (id+2) (id+3) ) are unaffected by what happened here (i.e. id+2 moving to time t)   
        // their collision times with their neighbors were absolute times and forces do not change
        
        /* update the collisions of (id+1) and (id+2) (or the boundary) */
        diagnose("right next coll");
        update_particle_task(id+1); 
        // there was a nontrivial collision with id+2 previously, which is now at particle id (due to swap) - delete it! 
        collisions.erase(particles[id].task);

        // if id-1 exists, update its collision with new id (which was id+1 before the swap) 
        if (id > 0) {
            /* get left neighbor particle to time t */
            update_particle(id-1, t);
            if (particles[id-1].x >= meanpos) {
                std::cout.precision(18);
                std::cout << "oops. " << id - 1 << " is at " << particles[id-1].x << " and meanpos is " << meanpos << std::endl;
                particles[id-1].x = meanpos - 0.00000001;
                std::cout << "now " << id - 1 << " is at " << particles[id-1].x << std::endl;
            }
            assert(particles[id-1].x < meanpos);

            /* and remove its collision */
            collisions.erase(particles[id-1].task);

            diagnose("left next coll");
        }
        else { /*id == 0 is a bit special - left boundary needs no time evolution */ 
            diagnose("left boudary update");
            collisions.erase(leftBoundaryTask);
        }
        /* update collision with particle id in any case */
        update_particle_task(id-1);

        // finally, and in any case, update necessarily nontrivial task for particle id 
        diagnose("update main particle");
        update_particle_task(id);
        diagnose("post");

    }

    //for (int i = 0; i < nParticles; ++i) std::cout << " " << particles[i].index<< " " << particles[i].t;
    //std::cout << std::endl<<"coll part was " << id << " at " << t;
    //std::cout << std::endl << std::endl;
}
