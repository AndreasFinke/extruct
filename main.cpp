#include "timer.h"
#include <random>
#include <tbb/tbb.h>
//#include <Eigen/Geometry>
#include <iostream>


//#include <pybind11/numpy.h>

#include "pcg32.h"


#include "main.h"
#include "multiverse.h"
#include "universe.h"
#include "background.h"


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


#include <memory>

#if PY == 0 

int main() {

    START_NAMED_TIMER("Big Bang")

    Background bg;
    bg.integrate();
    std::cout << bg.isIntegrated;
    Multiverse mv;
    PowerLaw * pl = new PowerLaw();

    for (int i = 0; i < 10; ++i) {
        mv.bang(1000, bg, pl);
    }

    STOP_TIMER

}
    //START_NAMED_TIMER("Evolution")

    //mv.evolveAll(0.5);
   
    //STOP_TIMER


    //Measurement * ps = new PowerSpectrumObs(0, 10, 1);
    //Measurement * d =  new DensityObs(100);


    ////mv.measure(0, d);
    //float * d_data = (float*)d->getResult();

    //std::cout << std::endl << std::endl << "Density = ";
    //for (int i = 0; i < 100; ++i) {
        ////std::cout << d_data[i] << " " << d_data[100+i] << std::endl;
    //}
    //std::cout <<std::endl << std::endl;

    //delete ps;

//}    

#endif

#if PY == 1 

PYBIND11_MODULE(extruct, m) {
    //py::class_<Universe>(m, "Universe")
        //.def(py::init<int>())
        //.def("draw", &Universe::draw)
        //.def("density", &Universe::density);
    py::class_<Multiverse>(m, "Multiverse")
        .def(py::init<>())
        .def("bang", &Multiverse::bang)
        .def("evolve", &Multiverse::evolve)
        .def("evolveAll", &Multiverse::evolveAll)
        .def("measure", &Multiverse::measure)
        .def("measureAll", &Multiverse::measureAll);
    py::class_<PowerSpectrum> powerspec(m, "PowerSpectrum");
    powerspec
        .def("eval", &PowerSpectrum::eval);
    py::class_<PowerLaw, PowerSpectrum>(m, "PowerLaw")
        .def(py::init<>())
        .def_readwrite("A", &PowerLaw::A);
    py::class_<Measurement> measurement(m, "Measurement");
    measurement
        .def("getResult", &Measurement::getResult)
        .def("measure", &Measurement::measure);
    py::class_<DensityObs, Measurement>(m, "DensityObs")
        .def(py::init<int>());
    py::class_<PhaseSpaceDensityObs, Measurement>(m, "PhaseSpaceDensityObs")
        .def(py::init<int>());
    py::class_<PowerSpectrumObs, Measurement>(m, "PowerSpectrumObs")
        .def(py::init<int, int, int>());
    py::class_<Background>(m, "Background")
        .def(py::init<>())
        .def(py::init<Float, Float, Float, Float, Float>())
        .def("integrate", &Background::integrate)
        .def("getScaleFactor", &Background::getScaleFactor)
        .def("getPhysTime", &Background::getPhysTime)
        .def("getD1", &Background::getD1)
        .def("getD2", &Background::getD2)
        .def("getD1d", &Background::getD1d);
        .def("getD2d", &Background::getD2d);
}

#endif
