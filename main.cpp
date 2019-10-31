#include "timer.h"
#include <random>
#include <tbb/tbb.h>
#include <iostream>

#include "ExactSum.h"
#include "sum.h"


//#include <pybind11/numpy.h>

#include "pcg32.h"


#include "main.h"
#include "multiverse.h"
#include "universe.h"
#include "background.h"

// there is still a sort of bug for less integration steps: periodic boundaries lead to D<0 in spline.h 
// understand why first branch never happens 
// connect points in phase space; find density. 
// power spectra as planned 
// replace ints by Longs
// do warnings
// coments as /* */ 


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

double fac(long n) {
    double ret = 1;
    for (Long i = 2; i <= n; ++i)
        ret *= i;
    return ret;
}

int main() {

    Float h = 0.7;
    Float omh2 = 0.12;
    Float obh2 = 0.0244;
    Float Om = (obh2 + omh2)/h/h;

    Float L = 20;
    
    //Background bg{100, 0, Om, h, 0.0};
    Background bg;
    bg.integrate();

    //std::cout <<  "Final tau = " << bg.getFinalTau() << std::endl;
    Multiverse mv(0);
    BBKS * pl = new BBKS(bg, L, true, false, obh2);
    pl->A = 1/(bg.getGrowth(0)*bg.getGrowth(0));



    DensityLin * d = new DensityLin(1024);

    //for (int i = 0; i < 10; ++i) {
        //mv.bang(1000, bg, pl, 1, 1);
    //}

    mv.bangAll(1024*8, bg, pl, L, 1, 1);
    //STOP_TIMER

    START_NAMED_TIMER("Evolution")

    mv.evolveAll(0);
    //std::cout << "first" << std::endl;
    //mv.evolve(9, 0);
    //std::cout << "second" << std::endl;
    //mv.evolve(8, 0);
    STOP_TIMER

    START_NAMED_TIMER("Density2")
    mv.measureAll(d);
    auto den = d->getResult();

   
    STOP_TIMER


    delete d;
    delete pl;
}

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
        .def(py::init<Long, bool>())
        .def("bang", &Multiverse::bang)
        .def("bangAll", &Multiverse::bangAll)
        .def("evolve", &Multiverse::evolve)
        .def("evolveAll", &Multiverse::evolveAll)
        .def("measure", &Multiverse::measure)
        .def("measureAll", &Multiverse::measureAll);
    py::class_<PowerSpectrum> powerspec(m, "PowerSpectrum");
    powerspec
        .def("eval_dimless", &PowerSpectrum::eval_dimless);
    py::class_<PowerLaw, PowerSpectrum>(m, "PowerLaw")
        .def(py::init<>())
        .def_readwrite("A", &PowerLaw::A);
    py::class_<BBKS, PowerSpectrum>(m, "BBKS")
        .def("eval", &BBKS::eval)
        .def(py::init<const Background&, Float, bool, bool, Float>())
        .def_readwrite("A", &BBKS::A)
        .def_readwrite("n", &BBKS::n)
        .def("set_dimless", &BBKS::set_dimless)
        .def("numModes", &BBKS::numModes); 
    py::class_<Measurement> measurement(m, "Measurement");
    measurement
        .def("getResult", &Measurement::getResult)
        .def("measure", &Measurement::measure);
    py::class_<DensityObs, Measurement>(m, "DensityObs")
        .def(py::init<int>());
    py::class_<DisplacementField, Measurement>(m, "DisplacementField")
        .def(py::init<int>());
    py::class_<DensityObs2, Measurement>(m, "DensityObs2")
        .def(py::init<int>());
    py::class_<DensityLin, Measurement>(m, "DensityLinObs")
        .def(py::init<int>());
    py::class_<DensitySPT<1>, Measurement>(m, "DensitySPT1")
        .def(py::init<int>());
    py::class_<DensitySPT<3>, Measurement>(m, "DensitySPT3")
        .def(py::init<int>());
    py::class_<DensitySPT<5>, Measurement>(m, "DensitySPT5")
        .def(py::init<int>());
    py::class_<DensitySPT<10>, Measurement>(m, "DensitySPT10")
        .def(py::init<int>());
    py::class_<DensitySPT<20>, Measurement>(m, "DensitySPT20")
        .def(py::init<int>());
    py::class_<DensitySPT<50>, Measurement>(m, "DensitySPT50")
        .def(py::init<int>());
    py::class_<DensitySPT<60>, Measurement>(m, "DensitySPT60")
        .def(py::init<int>());
    py::class_<DensitySPT<70>, Measurement>(m, "DensitySPT70")
        .def(py::init<int>());
    py::class_<DensitySPT<100>, Measurement>(m, "DensitySPT100")
        .def(py::init<int>());
    py::class_<DensitySPT<200>, Measurement>(m, "DensitySPT200")
        .def(py::init<int>());
    py::class_<DensitySPT<300>, Measurement>(m, "DensitySPT300")
        .def(py::init<int>());
    py::class_<DensitySPT<1000>, Measurement>(m, "DensitySPT1000")
        .def(py::init<int>());
    py::class_<CollisionObs, Measurement>(m, "CollisionObs")
        .def(py::init<int>());
    py::class_<PhaseSpaceDensityObs, Measurement>(m, "PhaseSpaceDensityObs")
        .def(py::init<int>());
    py::class_<PowerSpectrumObs, Measurement>(m, "PowerSpectrumObs")
        .def(py::init<int, int, Float, int>());
    py::class_<CorrelationFunctionObs, Measurement>(m, "CorrelationFunctionObs")
        .def(py::init<int, Float, bool>());
    py::class_<LinCorrelationFunctionObs, Measurement>(m, "LinCorrelationFunctionObs")
        .def(py::init<int, Float, bool>());
    py::class_<SPTCorrelationFunctionObs, Measurement>(m, "SPTCorrelationFunctionObs")
        .def(py::init<int, Float, bool>());
    py::class_<CorrelationFunction<DensitySPT<50>>, Measurement>(m, "SPT50CorrelationFunctionObs")
        .def(py::init<int, Float, bool>());
    py::class_<CorrelationFunction<DensitySPT<70>>, Measurement>(m, "SPT70CorrelationFunctionObs")
        .def(py::init<int, Float, bool>());
    py::class_<CorrelationFunction<DensitySPT<30>>, Measurement>(m, "SPT30CorrelationFunctionObs")
        .def(py::init<int, Float, bool>());
    py::class_<CorrelationFunction<DensitySPT<20>>, Measurement>(m, "SPT20CorrelationFunctionObs")
        .def(py::init<int, Float, bool>());
    py::class_<CorrelationFunction<DensitySPT<10>>, Measurement>(m, "SPT10CorrelationFunctionObs")
        .def(py::init<int, Float, bool>());
    py::class_<CorrelationFunction<DensitySPT<2>>, Measurement>(m, "SPT2CorrelationFunctionObs")
        .def(py::init<int, Float, bool>());
    py::class_<CorrelationFunction<DensitySPT<3>>, Measurement>(m, "SPT3CorrelationFunctionObs")
        .def(py::init<int, Float, bool>());
    py::class_<CorrelationFunction<DensitySPT<4>>, Measurement>(m, "SPT4CorrelationFunctionObs")
        .def(py::init<int, Float, bool>());
    //py::class_<PowerSpectrum3DObs, Measurement>(m, "PowerSpectrum3DObs")
        //.def(py::init<int, int>());
    py::class_<Background>(m, "Background")
        .def(py::init<>())
        .def(py::init<Float, Float, Float, Float, Float>())
        .def("integrate", &Background::integrate)
        .def("getScaleFactor", &Background::getScaleFactor)
        .def("getPhysTime", &Background::getPhysTime)
        .def("getD1", &Background::getD1)
        .def("getD2", &Background::getD2)
        .def("getD1d", &Background::getD1d)
        .def("getD2d", &Background::getD2d)
        .def("getPec", &Background::getPec)
        .def("getPecd", &Background::getPecd)
        .def("getGrowth", &Background::getGrowthOfZ);
}

#endif
