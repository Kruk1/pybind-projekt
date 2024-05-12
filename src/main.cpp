#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
#include <cmath>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)
#define _USE_MATH_DEFINES

int add(int i, int j) {
    return i + j;
}

float generate_signal(double frequency, double freq_samp, int samp)
{
    double faza = fmod(samp * 2 * M_PI * frequency / freq_samp, 2 * M_PI);
    return faza;
}

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: scikit_build_example

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";

    m.def("add", &add, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");

    m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");
    
    m.def("generate_signal", &generate_signal, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");
   
    m.attr("__version__") = "dev";
}
