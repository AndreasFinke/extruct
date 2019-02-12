#include "timer.h"
#include <random>
#include <tbb/tbb.h>
//#include <Eigen/Geometry>
#include <iostream>



//#include <pybind11/numpy.h>

#include "pcg32.h"


#include "main.h"
#include "universe.h"


struct Entry {
    Entry(Long idx, Float t) : t(t), idx(idx) {} 
    Float t;
    Long idx;
    bool operator<(const Entry& rhs) const {return t < rhs.t;}
    friend std::ostream& operator<<(std::ostream& out, const Entry& obj);
};

std::ostream& operator<<(std::ostream& out, const Entry& obj)
{
    return out << "entry_" << obj.t; 
}

template <template <typename, typename> class C, typename T, typename A>
void print_contents(const C<T, A>& v, const char * const separator = ", ")
{
    if(!v.empty())
    {
        std::copy(v.begin(),
                  --v.end(),
                  std::ostream_iterator<T>(std::cout, separator));
        std::cout << v.back() << "\n";
    }
}





int main() {

    START_NAMED_TIMER("Big Bang")
    Universe universe(0, 10);
    STOP_TIMER

    START_NAMED_TIMER("density")
    //universe.draw();
    //universe.draw();
    //universe.draw();
    int k = 0;
    for (; universe.latestTime < 400000; ++k) {
    //for (; k<200000; ++k) {
        double * den = universe.density();
        //std::cout << std::endl << std::endl;
        //for (int i = 0; i < universe.nParticles; ++i) {
            //std::cout << i << ": " << universe.get_particle_pos(i) << " (t=" << universe.get_particle_time(i) << ") ";
        //}

        //std::cout << std::endl;
        delete[] den; 
        universe.update_collision();
        //universe.synchronize();
        //universe.draw();
        //std::cout << std::endl << std::endl;
    }
    universe.synchronize();
    for (int i = 0; i < universe.nParticles; ++i) 
        std::cout << i << ": " << universe.get_particle_pos(i) << " (t=" << universe.get_particle_time(i) << ") ";

    std::cout << "integrated " << k << " collisions. " << std::endl;

    STOP_TIMER

    //Float a = pcg.nextFloat();

    //std::vector<Entry> q;
    //std::multiset<Entry> s;
    //const Long PARS = 400000;
    //const Long COLLS = 1000000;

    //for (Long i = 0; i < PARS; ++i) {
        //q.push_back(Entry(i, pcg.nextFloat()));
        //s.insert(Entry(i, pcg.nextFloat()));
    //}
    //std::sort(q.begin(), q.end());

    //START_NAMED_TIMER("vector")


    //for (Long j = 0; j < COLLS; ++j) {
        ////print_contents(q);
        //Float r = pcg.nextFloat();
        //q[0].t += r*r*r*0.01; 
        //auto pos = std::upper_bound(q.begin()+1, q.end(), q[0]);
        ////std::cout << pos-q.begin() << " ";
        ////auto test = q.begin();
        ////++test;
        ////std::cout << *test << " ";
        //std::rotate(q.begin(), q.begin()+1, pos);
        ////std::cout << *test << std::endl;
    //}


    //STOP_TIMER

    ////print_contents(q);

    //START_NAMED_TIMER("set")

    //for (Long j = 0; j < COLLS; ++j) {
        ////print_contents(q);
        //Float r = pcg.nextFloat();
        //Long idx = s.begin()->idx;
        //Float t = s.begin()->t;
        ////auto test = s.begin();
        ////++test;
        ////std::cout << *test << " ";
        //s.erase(s.begin());
        ////std::cout << *test << " ";
        //s.insert(Entry(idx, t+r*r*r*0.01));
        ////std::cout << *test << std::endl;
    //}
        
    //STOP_TIMER
}    

//PYBIND11_MODULE(extruct, m) {
    //py::class_<Universe>(m, "Universe")
        //.def(py::init<int>())
        //.def("draw", &Universe::draw)
        //.def("density", &Universe::density);
//}
//PYBIND11_MODULE(extruct, m) {
    //m.def("density", &density);
    //m.def("reset", &reset);
    //m.def("bigbang", &bigbang);
//}
